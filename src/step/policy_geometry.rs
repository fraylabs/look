//! Newtype wrappers that apply [`MeshingPolicy`] to STEP geometry without
//! modifying Truck.
//!
//! Truck's tessellation entry point takes a single scalar tolerance and
//! derives every curve's polyline and every surface's interior grid from it,
//! so a per-feature angular floor cannot be expressed through that scalar
//! alone. Rather than reimplement subdivision (which would have to re-derive
//! Truck's handling of partial arcs, full circles, reversed parameter ranges,
//! and periodic seam crossings), these wrappers forward *every* trait method
//! to the underlying STEP geometry unchanged — except
//! [`ParameterDivision1D::parameter_division`] /
//! [`ParameterDivision2D::parameter_division`], which compute an effective
//! tolerance (or, for surfaces, a minimum circumferential sample count) from
//! the policy and then delegate to the inner geometry. Truck's own machinery
//! does the actual subdivision, so the hard cases are handled exactly as
//! before.
//!
//! The wrappers preserve the inner geometry's identity, so the cylinder/cone
//! identification callbacks used by the formal recovery routes still see the
//! real `Surface`/`Curve3D` after a one-line unwrap. [`wrap_shell`] rebuilds
//! a compressed shell with wrapped edges and faces; because a compressed shell
//! samples each edge once by id, the shared-canonical-polyline invariant Truck
//! already relies on is kept.
//!
//! # Topology-safe eligibility
//!
//! Densifying a circular edge flips the constrained-Delaunay triangulator on
//! neighboring surfaces whose realization is density-sensitive — observed in
//! the NIST corpus as bidirectional sphere and B-spline face-recovery flips.
//! The cause is not the targeted cylinder/cone interiors but globally densifying
//! a *shared* circular boundary that also trims a sphere or a spline. So
//! "circular" alone is not sufficient eligibility.
//!
//! A circular edge is densified only when **every** incident face belongs to
//! the stable target set {plane, cylinder, cone}. A typical mechanical hole —
//! planar face ↔ circular edge ↔ cylindrical wall — stays eligible; a blended
//! cylinder ↔ sphere or cylinder ↔ spline transition keeps baseline sampling
//! on the shared edge. Cylindrical/conical interior floors are gated the same
//! way: a revolved face's circumferential floor applies only when its circular
//! boundaries are all eligible, avoiding a dense-interior/coarse-boundary
//! mismatch on mixed neighborhoods. Because the whole connected neighborhood
//! uses one policy per edge, adjacency stays crack-free.
//! A newly admitted shell containing a source-signed torus can include tori in
//! that target set only when the entire shell uses these analytic carriers.
//! Its positive companion tori share the same sampling policy; both native
//! angular axes receive compatible floors. Shells without signed tori and
//! mixed sphere/spline shells keep 0.38 torus sampling.

use std::collections::{HashMap, HashSet};

use truck_meshalgo::prelude::{
    BoundedCurve, BoundedSurface, D2, InnerSpace, ParameterDivision1D, ParameterDivision2D,
    ParameterRange, ParametricCurve, ParametricSurface, ParametricSurface3D, Point3, SPHint2D,
    SearchNearestParameter, SearchParameter, Vector3,
};
use truck_stepio::r#in::step_geometry::{Conic3D, Curve3D, ElementarySurface, Surface};
use truck_topology::compress::{CompressedEdge, CompressedFace, CompressedShell};

use crate::step::circular_arc::{decode_source_circle, decode_transformed_circle};
use crate::step::lattice::SplineAxisClosure;
use crate::step::meshing_policy::MeshingPolicy;

/// The stable target surface set: analytic surfaces whose constrained
/// triangulator is not density-sensitive. A circular edge is eligible for the
/// angular floor only when every incident face is one of these.
fn is_target_surface(surface: &Surface) -> bool {
    matches!(
        surface,
        Surface::ElementarySurface(ElementarySurface::Plane(_))
            | Surface::ElementarySurface(ElementarySurface::CylindricalSurface(_))
            | Surface::ElementarySurface(ElementarySurface::ConicalSurface(_))
    )
}

/// Whether a surface is a cylinder or cone (the revolved elementary surfaces
/// whose circumferential direction the policy can floor).
fn is_revolved_target(surface: &Surface) -> bool {
    matches!(
        surface,
        Surface::ElementarySurface(ElementarySurface::CylindricalSurface(_))
            | Surface::ElementarySurface(ElementarySurface::ConicalSurface(_))
    )
}

/// The certified radius of a STEP curve if it is a circle, else `None`.
///
/// A source-declared `circle` is decoded by family; a bare `ellipse` is
/// decoded by the exact Gram predicate, which refuses a genuinely non-circular
/// ellipse. Only certified circles can receive the angular floor.
fn circular_radius(curve: &Curve3D) -> Option<f64> {
    match curve {
        Curve3D::Conic(Conic3D::Circle(ellipse)) => decode_source_circle(ellipse)
            .ok()
            .map(|arc| arc.radius().get()),
        Curve3D::Conic(Conic3D::Ellipse(ellipse)) => decode_transformed_circle(ellipse)
            .ok()
            .map(|arc| arc.radius().get()),
        _ => None,
    }
}

/// A STEP edge curve carrying the meshing policy and a precomputed eligibility
/// flag.
///
/// Forwards [`ParametricCurve`] / [`BoundedCurve`] unchanged and overrides
/// [`ParameterDivision1D`] to apply the angular floor — but only when
/// [`eligible`](Self::new) is set, i.e. the edge is a certified circle whose
/// every incident face is in the target set. Ineligible edges (non-circles,
/// or circles shared with a sphere/spline/torus/etc.) divide at the baseline
/// linear tolerance, exactly as before.
#[derive(Clone, Debug)]
pub struct PolicyCurve {
    inner: Curve3D,
    policy: MeshingPolicy,
    eligible: bool,
}

impl PolicyCurve {
    /// Wrap a STEP curve with a policy and an eligibility flag.
    pub fn new(inner: Curve3D, policy: MeshingPolicy, eligible: bool) -> Self {
        Self {
            inner,
            policy,
            eligible,
        }
    }

    /// The underlying STEP curve, for the identification callbacks.
    pub fn inner(&self) -> &Curve3D {
        &self.inner
    }
}

impl ParametricCurve for PolicyCurve {
    type Point = Point3;
    type Vector = Vector3;
    #[inline]
    fn subs(&self, t: f64) -> Self::Point {
        self.inner.subs(t)
    }
    #[inline]
    fn der(&self, t: f64) -> Self::Vector {
        self.inner.der(t)
    }
    #[inline]
    fn der2(&self, t: f64) -> Self::Vector {
        self.inner.der2(t)
    }
    #[inline]
    fn der_n(&self, n: usize, t: f64) -> Self::Vector {
        self.inner.der_n(n, t)
    }
    #[inline]
    fn parameter_range(&self) -> ParameterRange {
        self.inner.parameter_range()
    }
    #[inline]
    fn try_range_tuple(&self) -> Option<(f64, f64)> {
        self.inner.try_range_tuple()
    }
    #[inline]
    fn period(&self) -> Option<f64> {
        self.inner.period()
    }
}

impl BoundedCurve for PolicyCurve {
    #[inline]
    fn evaluation_range(&self) -> (f64, f64) {
        self.inner.evaluation_range()
    }

    #[inline]
    fn basis_is_partition_of_unity(&self, t: f64) -> bool {
        self.inner.basis_is_partition_of_unity(t)
    }
}

impl ParameterDivision1D for PolicyCurve {
    type Point = Point3;
    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
        let eff = if self.eligible {
            match circular_radius(&self.inner) {
                Some(radius) => self.policy.effective_curve_tolerance(tol, radius),
                None => tol,
            }
        } else {
            // Ineligible: the linear tolerance is the only mechanism, and it is
            // already world-space. Forward unchanged.
            tol
        };
        self.inner.parameter_division(range, eff)
    }
}

/// A STEP surface carrying the meshing policy and a precomputed
/// interior-floor eligibility flag.
///
/// Forwards [`ParametricSurface`] / [`SearchParameter`] / [`SearchNearestParameter`]
/// unchanged and overrides [`ParameterDivision2D`] to hold the circumferential
/// direction of cylindrical and conical surfaces to the policy's angular floor
/// — but only when [`interior_eligible`](Self::new) is set, i.e. the face is a
/// cylinder/cone or torus in a recovered analytic shell whose circular boundaries are eligible. Otherwise the
/// surface divides at the baseline tolerance, preserving the
/// dense-interior/coarse-boundary match.
/// A per-axis quotient map from cover (deck-lifted) coordinates to the native
/// evaluator interval of a source-certified closed spline axis.
///
/// The topology reasons about unwrapped cover UV; the generic spline evaluator
/// supports only its native `[a, b]`. On a certified axis the physical
/// representative is
///
/// ```text
/// x_native = a + (x_cover - a).rem_euclid(P),   P = b - a
/// ```
///
/// which is continuous across the seam because the certified axis satisfies the
/// position and first-derivative seam identification. Non-periodic and analytic
/// axes never carry a quotient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuotientAxis {
    /// The native evaluation interval `[a, b]` on this axis.
    pub interval: (f64, f64),
    /// The period `P = b - a`.
    pub period: f64,
}

impl QuotientAxis {
    /// The native evaluator representative of a cover coordinate on this axis.
    #[inline]
    pub fn map(&self, cover: f64) -> f64 {
        let (a, _) = self.interval;
        a + (cover - a).rem_euclid(self.period)
    }
}

/// A STEP surface carrying the meshing policy and a precomputed
/// interior-floor eligibility flag.
///
/// Forwards [`ParametricSurface`] / [`SearchParameter`] / [`SearchNearestParameter`]
/// unchanged and overrides [`ParameterDivision2D`] to hold the circumferential
/// direction of cylindrical and conical surfaces to the policy's angular floor
/// — but only when [`interior_eligible`](Self::new) is set, i.e. the face is a
/// cylinder/cone or torus in a recovered analytic shell whose circular boundaries are eligible. Otherwise the
/// surface divides at the baseline tolerance, preserving the
/// dense-interior/coarse-boundary match.
///
/// # Cover UV vs native evaluator UV
///
/// When the composition layer attaches a source-certified closed spline axis
/// ([`with_closure`](Self::with_closure)), the physical evaluator methods
/// ([`subs`](Self::subs) and the derivative/normal family) map the incoming
/// cover UV to the native evaluator representative through the per-axis
/// quotient before forwarding. Topology-facing methods (`parameter_range`,
/// `try_range_tuple`, `u_period`, `v_period`) and the UV the tessellator stores
/// are never rewritten: the quotient is a *selection of the evaluator
/// representative*, never a topology identification. Analytic periodic surfaces
/// carry no quotient and forward byte-identically.
#[derive(Clone, Debug)]
pub struct PolicySurface {
    inner: Surface,
    policy: MeshingPolicy,
    interior_eligible: bool,
    /// The source-declared spline-axis closure attached by the composition
    /// layer, when the STEP table names this face's support surface. Carried
    /// as already-established provenance: no STEP inference happens here.
    source_closure: Option<SplineAxisClosure>,
    /// The cover→native evaluator quotient on the `u` axis, when certified.
    u_quotient: Option<QuotientAxis>,
    /// The cover→native evaluator quotient on the `v` axis, when certified.
    v_quotient: Option<QuotientAxis>,
    native_inverse_recovery: bool,
    native_control_hull: std::sync::OnceLock<Option<(Point3, Point3)>>,
}

impl PolicySurface {
    /// Wrap a STEP surface with a policy and an interior-floor eligibility flag.
    pub fn new(inner: Surface, policy: MeshingPolicy, interior_eligible: bool) -> Self {
        Self {
            inner,
            policy,
            interior_eligible,
            source_closure: None,
            u_quotient: None,
            v_quotient: None,
            native_inverse_recovery: false,
            native_control_hull: std::sync::OnceLock::new(),
        }
    }

    /// Wrap a STEP surface with a policy, an interior-floor eligibility flag,
    /// and the source-declared spline-axis closure of its support surface.
    ///
    /// The per-axis quotient descriptors are derived once from the same lattice
    /// and `evaluation_range` the certification theorem used, so the evaluator
    /// map carries already-established provenance rather than re-deriving
    /// periodicity numerically.
    pub fn with_closure(
        inner: Surface,
        policy: MeshingPolicy,
        interior_eligible: bool,
        source_closure: Option<SplineAxisClosure>,
    ) -> Self {
        let (u_interval, v_interval) =
            crate::step::lattice::spline_quotient_axes(&inner, source_closure);
        let u_quotient = u_interval.map(|(a, b)| QuotientAxis {
            interval: (a, b),
            period: b - a,
        });
        let v_quotient = v_interval.map(|(a, b)| QuotientAxis {
            interval: (a, b),
            period: b - a,
        });
        Self {
            inner,
            policy,
            interior_eligible,
            source_closure,
            u_quotient,
            v_quotient,
            native_inverse_recovery: false,
            native_control_hull: std::sync::OnceLock::new(),
        }
    }

    /// A search-only retry for an ordinary spline whose original trim refused.
    /// The carrier, lattice and physical boundary samples stay unchanged.
    pub fn native_inverse_retry(&self) -> Option<Self> {
        if self.u_quotient.is_some()
            || self.v_quotient.is_some()
            || !matches!(
                self.inner,
                Surface::BSplineSurface(_) | Surface::NurbsSurface(_)
            )
        {
            return None;
        }
        let mut retry = self.clone();
        retry.native_inverse_recovery = true;
        Some(retry)
    }

    /// The underlying STEP surface, for the identification callbacks.
    pub fn inner(&self) -> &Surface {
        &self.inner
    }

    /// The source-declared spline-axis closure, as attached by the composition
    /// layer. `None` when no STEP table entry names this face's support surface.
    pub fn source_closure(&self) -> Option<SplineAxisClosure> {
        self.source_closure
    }

    /// The certified cover→native evaluator quotient on the `u` axis, if any.
    pub fn u_quotient(&self) -> Option<QuotientAxis> {
        self.u_quotient
    }

    /// The certified cover→native evaluator quotient on the `v` axis, if any.
    pub fn v_quotient(&self) -> Option<QuotientAxis> {
        self.v_quotient
    }

    /// The native evaluator representative of a cover UV pair.
    ///
    /// Non-periodic axes pass through unchanged; analytic periodic surfaces
    /// carry no quotient, so they too pass through unchanged (their evaluators
    /// are globally periodic and must not be routed through a mod map).
    #[inline]
    fn native_uv(&self, u: f64, v: f64) -> (f64, f64) {
        (
            self.u_quotient.map_or(u, |q| q.map(u)),
            self.v_quotient.map_or(v, |q| q.map(v)),
        )
    }

    fn inverse_is_in_native_domain(&self, uv: (f64, f64)) -> bool {
        let ((u0, u1), (v0, v1)) = match &self.inner {
            Surface::BSplineSurface(s) => BoundedSurface::evaluation_range(s),
            Surface::NurbsSurface(s) => BoundedSurface::evaluation_range(s),
            _ => return true,
        };
        let in_range = |x: f64, a: f64, b: f64| {
            let margin = (b - a).abs().max(1.0) * 1.0e-6;
            x >= a - margin && x <= b + margin
        };
        (self.u_quotient.is_some() || in_range(uv.0, u0, u1))
            && (self.v_quotient.is_some() || in_range(uv.1, v0, v1))
    }

    fn native_control_hull_distance(&self, point: Point3) -> Option<f64> {
        let (min, max) = (*self
            .native_control_hull
            .get_or_init(|| self.compute_native_control_hull()))?;
        Some(
            (0..3)
                .map(|axis| {
                    let delta = (min[axis] - point[axis])
                        .max(point[axis] - max[axis])
                        .max(0.0);
                    delta * delta
                })
                .sum::<f64>()
                .sqrt(),
        )
    }

    fn compute_native_control_hull(&self) -> Option<(Point3, Point3)> {
        let points: Vec<Point3> = match &self.inner {
            Surface::BSplineSurface(s) => s.control_points().iter().flatten().copied().collect(),
            Surface::NurbsSurface(s) => {
                let controls = s.control_points().iter().flatten();
                if controls.clone().any(|p| !p.w.is_finite() || p.w <= 0.0) {
                    return None;
                }
                controls
                    .map(|p| Point3::new(p.x / p.w, p.y / p.w, p.z / p.w))
                    .collect()
            }
            _ => return None,
        };
        if points.is_empty()
            || points
                .iter()
                .any(|p| !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite())
        {
            return None;
        }
        let mut min = points[0];
        let mut max = min;
        for p in points {
            for axis in 0..3 {
                min[axis] = min[axis].min(p[axis]);
                max[axis] = max[axis].max(p[axis]);
            }
        }
        Some((min, max))
    }

    fn inverse_is_strictly_native(&self, uv: (f64, f64)) -> bool {
        let ((u0, u1), (v0, v1)) = match &self.inner {
            Surface::BSplineSurface(s) => BoundedSurface::evaluation_range(s),
            Surface::NurbsSurface(s) => BoundedSurface::evaluation_range(s),
            _ => return true,
        };
        uv.0 >= u0 && uv.0 <= u1 && uv.1 >= v0 && uv.1 <= v1
    }

    fn inverse_has_regular_jacobian(&self, uv: (f64, f64), spans: (f64, f64)) -> bool {
        // Normalize derivatives to the native domain, so the rank decision is
        // invariant to carrier units and independent parameter rescaling.
        // A nearly collapsed edge has many indistinguishable inverses; choosing
        // another root there can change trim order without moving its point.
        let du = self.inner.uder(uv.0, uv.1) * spans.0;
        let dv = self.inner.vder(uv.0, uv.1) * spans.1;
        let scale = du.magnitude2() + dv.magnitude2();
        let area2 = du.cross(dv).magnitude2();
        scale.is_finite() && area2.is_finite() && area2 > 64.0 * f64::EPSILON * scale * scale
    }

    fn prefer_native_inverse(&self, point: Point3, uv: (f64, f64), trials: usize) -> (f64, f64) {
        if self.u_quotient.is_some()
            || self.v_quotient.is_some()
            || self.inverse_is_in_native_domain(uv)
        {
            return uv;
        }
        // An exterior Newton root can realize the boundary point while its
        // extrapolated chart produces unrelated interior geometry. Prefer an
        // equally accurate native-domain root, without rejecting legacy faces
        // for which the source admits only an exterior inverse.
        let residual = (self.inner.subs(uv.0, uv.1) - point).magnitude();
        // On the active knot domain, polynomial splines and positive-weight
        // rational splines lie in their control hull. A point farther from
        // that hull than the accepted residual cannot have a qualifying root.
        // New roots must be strictly native for this bound to apply; the
        // existing padded-domain acceptance of the original root is unchanged.
        if self
            .native_control_hull_distance(point)
            .is_some_and(|distance| distance > residual + truck_meshalgo::prelude::TOLERANCE)
        {
            return uv;
        }
        let ((u0, u1), (v0, v1)) = match &self.inner {
            Surface::BSplineSurface(s) => BoundedSurface::evaluation_range(s),
            Surface::NurbsSurface(s) => BoundedSurface::evaluation_range(s),
            _ => return uv,
        };
        if ![u0, u1, v0, v1].iter().all(|x| x.is_finite()) || u1 <= u0 || v1 <= v0 {
            return uv;
        }
        let spans = (u1 - u0, v1 - v0);
        if !self.inverse_has_regular_jacobian(uv, spans) {
            return uv;
        }
        // Keep the original inverse as a locality hint, but start Newton on
        // the native chart. Only the starting guess is projected: the result
        // must independently satisfy the domain and physical residual checks.
        // A central start is a cheap second attempt. These two starts need not
        // reach another turn of a long spline strip, so failure must still
        // admit the native presearch and the source knot-span cells below.
        let seeds = [
            (uv.0.clamp(u0, u1), uv.1.clamp(v0, v1)),
            ((u0 + u1) / 2.0, (v0 + v1) / 2.0),
        ];
        let original_normal =
            (self.inner.uder(uv.0, uv.1) * spans.0).cross(self.inner.vder(uv.0, uv.1) * spans.1);
        let try_hint = |hint, preserve_orientation: bool| {
            let candidate = self
                .inner
                .search_parameter(point, hint, trials)
                .or_else(|| self.inner.search_nearest_parameter(point, hint, trials));
            if let Some(candidate) = candidate {
                if self.inverse_is_strictly_native(candidate)
                    && (self.inner.subs(candidate.0, candidate.1) - point).magnitude()
                        <= residual + truck_meshalgo::prelude::TOLERANCE
                    && self.inverse_has_regular_jacobian(candidate, spans)
                    && (!preserve_orientation
                        || original_normal.dot(
                            (self.inner.uder(candidate.0, candidate.1) * spans.0)
                                .cross(self.inner.vder(candidate.0, candidate.1) * spans.1),
                        ) > 0.0)
                {
                    return Some(candidate);
                }
            }
            None
        };
        if let Some(candidate) = seeds
            .into_iter()
            .find_map(|seed| try_hint(Some(seed), false))
        {
            return candidate;
        }
        // A successful trim owns its existing chart. Remote roots become
        // eligible only on the kernel's typed, failed-face retry path.
        if !self.native_inverse_recovery {
            return uv;
        }
        // Build the broader starts only after the local attempts fail. Keep
        // the same residual, native-domain and Jacobian requirements: the
        // fallback finds a chart representative, never moves a source point.
        // Unlike the established local starts, a remote cell can lie on a
        // different regular sheet of a folded carrier. Equal positions and
        // full-rank derivatives do not authorize reversing that sheet's
        // orientation; preserve the original normal hemisphere here.
        std::iter::once(None)
            .chain(self.inner.search_parameter_seeds().into_iter().map(Some))
            .find_map(|hint| try_hint(hint, true))
            .unwrap_or(uv)
    }

    /// Whether a final parameter-inverse result is a valid representative.
    ///
    /// The strict final-range rule exists to keep a *source-certified periodic
    /// policy surface* epistemically coherent: its B-spline inverse can converge
    /// to a spurious stationary root outside the native domain (Newton escapes
    /// past the knot ends), and the certified axis carries deck-lifted topology
    /// while the ordinary companion axis stays bounded. On a certified-closed
    /// axis the deck-equivalent (cover-lifted) values are legal — the quotient
    /// normalizes them during evaluation — so no bound applies there; the
    /// ordinary companion axis must stay within its `evaluation_range`, up to a
    /// justified relative tolerance.
    ///
    /// Ordinary B-spline/NURBS surfaces with **no certified periodic quotient**
    /// preserve legacy inverse semantics: their tessellation has always
    /// accepted the inverse result as-is, and imposing the range rule on them
    /// turns previously-rendering faces into `BoundaryProjectionFailed`.
    /// Analytic surfaces are globally periodic or unbounded and always accept.
    pub fn accept_inverse_result(&self, uv: (f64, f64)) -> bool {
        // The strict range rule is a periodic-policy instrument. A face with no
        // certified periodic axis on either axis is an ordinary spline and keeps
        // its legacy inverse behavior (NIST #1167-wave regression witness:
        // 59 ordinary B-spline/NURBS faces lost to `BoundaryProjectionFailed`
        // because the guard rejected their pre-existing inverse results).
        if self.u_quotient.is_none() && self.v_quotient.is_none() {
            return true;
        }
        let ((u0, u1), (v0, v1)) = match &self.inner {
            Surface::BSplineSurface(spline) => BoundedSurface::evaluation_range(spline),
            Surface::NurbsSurface(spline) => BoundedSurface::evaluation_range(spline),
            _ => return true,
        };
        let tol = |(a, b): (f64, f64)| (b - a).abs().max(1.0) * 1.0e-6;
        let u_ok = self.u_quotient.is_some() || {
            let t = tol((u0, u1));
            uv.0 >= u0 - t && uv.0 <= u1 + t
        };
        let v_ok = self.v_quotient.is_some() || {
            let t = tol((v0, v1));
            uv.1 >= v0 - t && uv.1 <= v1 + t
        };
        u_ok && v_ok
    }
}

impl ParametricSurface for PolicySurface {
    type Point = Point3;
    type Vector = Vector3;
    #[inline]
    fn subs(&self, u: f64, v: f64) -> Self::Point {
        let (u, v) = self.native_uv(u, v);
        self.inner.subs(u, v)
    }
    #[inline]
    fn uder(&self, u: f64, v: f64) -> Self::Vector {
        let (u, v) = self.native_uv(u, v);
        self.inner.uder(u, v)
    }
    #[inline]
    fn vder(&self, u: f64, v: f64) -> Self::Vector {
        let (u, v) = self.native_uv(u, v);
        self.inner.vder(u, v)
    }
    #[inline]
    fn uuder(&self, u: f64, v: f64) -> Self::Vector {
        let (u, v) = self.native_uv(u, v);
        self.inner.uuder(u, v)
    }
    #[inline]
    fn uvder(&self, u: f64, v: f64) -> Self::Vector {
        let (u, v) = self.native_uv(u, v);
        self.inner.uvder(u, v)
    }
    #[inline]
    fn vvder(&self, u: f64, v: f64) -> Self::Vector {
        let (u, v) = self.native_uv(u, v);
        self.inner.vvder(u, v)
    }
    #[inline]
    fn der_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Self::Vector {
        let (u, v) = self.native_uv(u, v);
        self.inner.der_mn(m, n, u, v)
    }
    #[inline]
    fn parameter_range(&self) -> (ParameterRange, ParameterRange) {
        self.inner.parameter_range()
    }
    #[inline]
    fn try_range_tuple(&self) -> (Option<(f64, f64)>, Option<(f64, f64)>) {
        self.inner.try_range_tuple()
    }
    #[inline]
    fn u_period(&self) -> Option<f64> {
        self.inner.u_period()
    }
    #[inline]
    fn v_period(&self) -> Option<f64> {
        self.inner.v_period()
    }
}

impl ParametricSurface3D for PolicySurface {
    // Forward the normal-family methods rather than inheriting the trait's
    // `uder × vder` defaults: the inner `Surface`'s derived `ParametricSurface3D`
    // delegates to per-variant impls, and a variant may carry a normal-family
    // override. Using the defaults here would diverge from the unwrapped path
    // and perturb borderline CDTs. The native-UV quotient applies here too.
    #[inline]
    fn normal(&self, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.native_uv(u, v);
        self.inner.normal(u, v)
    }
    #[inline]
    fn normal_uder(&self, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.native_uv(u, v);
        self.inner.normal_uder(u, v)
    }
    #[inline]
    fn normal_vder(&self, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.native_uv(u, v);
        self.inner.normal_vder(u, v)
    }
}

impl ParameterDivision2D for PolicySurface {
    fn parameter_division(
        &self,
        range: ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        // A source-certified periodic spline axis must be subdivided in one
        // locally coherent physical chart — its native fundamental domain — or
        // the recursive bilinear error estimator reads seam-adjacent corners as
        // numerically distant and refines every spanning cell to
        // `MAX_DIVISION_CELLS`. Subdivide the coherent chart, then map the
        // division back onto the unwrapped cover interval. Non-certified
        // surfaces keep the existing policy floor path unchanged.
        if self.u_quotient.is_some() || self.v_quotient.is_some() {
            // Counterfactual gate: `LOOK_NO_QUOTIENT_DIVISION` disables Stage F
            // only (the cover range is subdivided raw), for one-stage
            // attribution of a census regression.
            if std::env::var_os("LOOK_NO_QUOTIENT_DIVISION").is_none() {
                return self.quotient_parameter_division(range, tol);
            }
        }
        if !self.interior_eligible {
            return self.inner.parameter_division(range, tol);
        }
        // `RevolutedCurve::subs(u, v)` rotates by `v`, so `v` is the
        // circumferential (revolution) direction and `u` is the axial (generatrix)
        // direction. The linear tolerance, capped by the absolute deflection,
        // governs the axial direction and provides the base circumferential grid;
        // the angular floor then lifts the circumferential count if the linear
        // term left it below the policy's proportional share of a revolution.
        let capped_tol = tol.min(self.policy.maximum_absolute_deflection);
        let (mut udiv, mut vdiv) = self.inner.parameter_division(range, capped_tol);
        if matches!(
            self.inner,
            Surface::ElementarySurface(ElementarySurface::ToroidalSurface(_))
        ) {
            // Only tori in recovered analytic shells can be eligible here. Both native angular
            // axes must follow their refined shared circular boundaries.
            for (division, (min, max)) in [(&mut udiv, range.0), (&mut vdiv, range.1)] {
                let needed = self.policy.surface_floor_segments(max - min);
                if needed >= 2 && division.len() < needed {
                    *division = uniform_linspace(min, max, needed);
                }
            }
            return (udiv, vdiv);
        }
        let (v_min, v_max) = range.1;
        let needed = self.policy.surface_floor_segments(v_max - v_min);
        if needed >= 2 && vdiv.len() < needed {
            (udiv, uniform_linspace(v_min, v_max, needed))
        } else {
            (udiv, vdiv)
        }
    }
}

impl PolicySurface {
    /// Subdivide a cover interval whose certified periodic axes are folded to
    /// their native fundamental domain, then map the division back to the
    /// unwrapped cover interval.
    ///
    /// The physical content of a cover interval is period-invariant: translating
    /// the whole interval by a deck step changes only which copy each native
    /// sample lands in, never how many samples the physical patch needs. Mapping
    /// through the native chart is what makes tessellation complexity independent
    /// of deck offset.
    fn quotient_parameter_division(
        &self,
        range: ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let ((u0, u1), (v0, v1)) = range;
        // Certified axes are subdivided over their native fundamental domain
        // `[a, b]` (one coherent chart, no seam inside); the other axes keep the
        // requested cover interval.
        let (u_lo, u_hi) = self.u_quotient.map_or((u0, u1), |q| q.interval);
        let (v_lo, v_hi) = self.v_quotient.map_or((v0, v1), |q| q.interval);
        let (udiv_native, vdiv_native) = self
            .inner
            .parameter_division(((u_lo, u_hi), (v_lo, v_hi)), tol);
        let udiv = match self.u_quotient {
            Some(q) => map_division_to_cover((u0, u1), q.interval.0, q.period, &udiv_native),
            None => udiv_native,
        };
        let vdiv = match self.v_quotient {
            Some(q) => map_division_to_cover((v0, v1), q.interval.0, q.period, &vdiv_native),
            None => vdiv_native,
        };
        (udiv, vdiv)
    }
}

/// Map a subdivision of the native fundamental domain `[base, base + period]`
/// back onto an unwrapped cover interval, replicating each native sample at
/// every deck copy whose window intersects the cover interval and dropping the
/// copies that fall outside it.
///
/// The returned sequence is ascending, contains both cover endpoints, and its
/// length is invariant to translating the cover interval by whole periods —
/// the property the quotient subdivision theorem requires.
fn map_division_to_cover(
    cover: (f64, f64),
    base: f64,
    period: f64,
    native_div: &[f64],
) -> Vec<f64> {
    let (c0, c1) = cover;
    if c1 <= c0 {
        return native_div.to_vec();
    }
    let k_lo = ((c0 - base) / period).floor() as i64;
    let k_hi = ((c1 - base) / period).floor() as i64;
    let mut out: Vec<f64> = Vec::new();
    for k in k_lo..=k_hi {
        let offset = k as f64 * period;
        for &x in native_div {
            let cover_x = x + offset;
            if cover_x >= c0 && cover_x <= c1 {
                out.push(cover_x);
            }
        }
    }
    out.push(c0);
    out.push(c1);
    out.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    out.dedup_by(|a, b| (*a - *b).abs() < 1.0e-9);
    out
}

impl SearchParameter<D2> for PolicySurface {
    type Point = Point3;
    #[inline]
    fn search_parameter<H: Into<SPHint2D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        // On a source-certified periodic policy surface the inverse can converge
        // to a spurious root outside the native domain of the ordinary companion
        // axis (e.g. `u = 3.63` for a true `u = 0.765` on `[0,1]`). Rejecting the
        // final result here lets the existing fallback chain (hintless /
        // structural seeds) find the in-domain root; it never clamps an exterior
        // value onto the boundary. For an ordinary spline, prefer an equally
        // accurate native root before retaining the legacy result
        // (`accept_inverse_result` applies no bound).
        let uv = self.inner.search_parameter(point, hint, trials)?;
        let uv = self.prefer_native_inverse(point, uv, trials);
        self.accept_inverse_result(uv).then_some(uv)
    }
    // Forward the structure-derived seeds (e.g. a B-spline's knot-span starts)
    // rather than inheriting the empty default. The default drops the piecewise
    // geometry's own inverse hints, starving the projection Newton iteration
    // and turning faces that render on the unwrapped path into
    // `NoSurfaceProduced`.
    #[inline]
    fn search_parameter_seeds(&self) -> Vec<(f64, f64)> {
        self.inner.search_parameter_seeds()
    }
}

impl SearchNearestParameter<D2> for PolicySurface {
    type Point = Point3;
    #[inline]
    fn search_nearest_parameter<H: Into<SPHint2D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let uv = self.inner.search_nearest_parameter(point, hint, trials)?;
        let uv = self.prefer_native_inverse(point, uv, trials);
        self.accept_inverse_result(uv).then_some(uv)
    }
}

/// `n` parameters evenly spaced over `[a, b]`, endpoints inclusive. `n >= 2`.
fn uniform_linspace(a: f64, b: f64, n: usize) -> Vec<f64> {
    debug_assert!(n >= 2);
    let step = (b - a) / (n - 1) as f64;
    (0..n).map(|i| a + step * i as f64).collect()
}

/// Per-edge and per-face eligibility for the topology-safe angular floor.
///
/// [`edge_eligible`](Self::edge_eligible) is true for a canonical topological
/// edge exactly when it is a certified circle and every face incident on it is
/// in the target set {plane, cylinder, cone}. [`face_interior_eligible`]
/// `(Self::face_interior_eligible)` is true for a cylinder/cone face exactly
/// when every circular edge in its boundary is edge-eligible, so the
/// circumferential interior floor never creates a dense-interior/coarse-boundary
/// mismatch on a mixed neighborhood.
struct Eligibility {
    edge_eligible: Vec<bool>,
    face_interior_eligible: Vec<bool>,
}

impl Eligibility {
    #[cfg(test)]
    fn compute(shell: &CompressedShell<Point3, Curve3D, Surface>) -> Self {
        Self::with_signed_tori(shell, &HashSet::new())
    }

    fn with_signed_tori(
        shell: &CompressedShell<Point3, Curve3D, Surface>,
        signed_torus_surfaces: &HashSet<u64>,
    ) -> Self {
        // A signed torus makes this shell newly admissible. Refine its whole
        // analytic torus neighborhood, including positive companion carriers;
        // refining only the signed face leaves their shared outer outline at
        // the coarse assembly tolerance. Mixed spline/sphere shells keep the
        // old sampling because their triangulation is density-sensitive.
        let recovered_analytic_shell = shell.faces.iter().any(|face| {
            face.provenance
                .surface_id
                .is_some_and(|id| signed_torus_surfaces.contains(&id.get()))
                && matches!(
                    face.surface,
                    Surface::ElementarySurface(ElementarySurface::ToroidalSurface(_))
                )
        }) && shell.faces.iter().all(|face| {
            is_target_surface(&face.surface)
                || matches!(
                    face.surface,
                    Surface::ElementarySurface(ElementarySurface::ToroidalSurface(_))
                )
        });
        let recovered_torus_faces: Vec<bool> = shell
            .faces
            .iter()
            .map(|face| {
                recovered_analytic_shell
                    && matches!(
                        face.surface,
                        Surface::ElementarySurface(ElementarySurface::ToroidalSurface(_))
                    )
            })
            .collect();
        // Map each edge index to the faces that reference it. A face references
        // an edge once per boundary loop occurrence; collecting duplicates is
        // harmless for the all-targets check.
        let mut edge_faces: Vec<Vec<usize>> = vec![Vec::new(); shell.edges.len()];
        for (face_index, face) in shell.faces.iter().enumerate() {
            for boundary in &face.boundaries {
                for edge_index in boundary {
                    edge_faces[edge_index.index].push(face_index);
                }
            }
        }

        let edge_eligible: Vec<bool> = shell
            .edges
            .iter()
            .enumerate()
            .map(|(e, edge)| {
                // Must be a certified circle. `eligible` carrying a non-circle
                // would be harmless (it forwards `tol` unchanged) but computing
                // it here keeps the invariant that eligible ⟹ circle explicit.
                if circular_radius(&edge.curve).is_none() {
                    return false;
                }
                edge_faces[e].iter().all(|&f| {
                    is_target_surface(&shell.faces[f].surface) || recovered_torus_faces[f]
                })
            })
            .collect();

        let face_interior_eligible: Vec<bool> = shell
            .faces
            .iter()
            .enumerate()
            .map(|(f, face)| {
                if !is_revolved_target(&face.surface) && !recovered_torus_faces[f] {
                    return false;
                }
                // The interior floor applies only when every circular boundary
                // edge of this face is edge-eligible. A circular edge that is
                // shared with a sphere/spline/etc. is ineligible, so the face
                // keeps its baseline interior division to stay matched with
                // that edge's baseline boundary sampling.
                let mut all_circular_eligible = true;
                for boundary in &face.boundaries {
                    for edge_index in boundary {
                        let e = edge_index.index;
                        if circular_radius(&shell.edges[e].curve).is_some() && !edge_eligible[e] {
                            all_circular_eligible = false;
                        }
                    }
                }
                all_circular_eligible
            })
            .collect();

        Self {
            edge_eligible,
            face_interior_eligible,
        }
    }
}

/// Rebuild a compressed shell with policy-wrapped edge curves and surfaces.
///
/// Vertices, topology (edge indices, boundaries, orientation), and face
/// provenance are carried through unchanged; only the curve and surface
/// geometries are wrapped. Eligibility is computed once from the shell's
/// topology and carried into each wrapper, so a shared edge has a single
/// density for all incident faces (crack-free) and the angular floor applies
/// only inside topology-safe plane/cylinder/cone neighborhoods.
pub fn wrap_shell(
    shell: CompressedShell<Point3, Curve3D, Surface>,
    policy: MeshingPolicy,
) -> CompressedShell<Point3, PolicyCurve, PolicySurface> {
    wrap_shell_with_closure(shell, policy, &HashMap::new())
}

/// As [`wrap_shell`], additionally attaching each face's source-declared
/// spline-axis closure, keyed by the STEP face definition or surface entity
/// id its provenance names. Face entries account for converted parameter axes.
///
/// The closure map is built by the composition layer from the STEP table
/// ([`crate::step::lattice::spline_closure_map`]); this wrapper only carries the
/// established metadata into each [`PolicySurface`]. A face whose support
/// surface is not a spline, or whose provenance names no surface, receives
/// `None`.
pub fn wrap_shell_with_closure(
    shell: CompressedShell<Point3, Curve3D, Surface>,
    policy: MeshingPolicy,
    closure_map: &HashMap<u64, SplineAxisClosure>,
) -> CompressedShell<Point3, PolicyCurve, PolicySurface> {
    wrap_shell_with_source_metadata(shell, policy, closure_map, &HashSet::new())
}

/// Carry source closure and signed-major torus identity into the shared policy.
/// Newly admitted, entirely analytic shells containing a signed torus use the
/// same angular floor on their torus carriers and stable neighboring faces.
/// Other torus neighborhoods keep 0.38 sampling.
pub fn wrap_shell_with_source_metadata(
    shell: CompressedShell<Point3, Curve3D, Surface>,
    policy: MeshingPolicy,
    closure_map: &HashMap<u64, SplineAxisClosure>,
    signed_torus_surfaces: &HashSet<u64>,
) -> CompressedShell<Point3, PolicyCurve, PolicySurface> {
    let eligibility = Eligibility::with_signed_tori(&shell, signed_torus_surfaces);
    CompressedShell {
        vertices: shell.vertices,
        edges: shell
            .edges
            .into_iter()
            .enumerate()
            .map(|(e, edge)| CompressedEdge {
                vertices: edge.vertices,
                curve: PolicyCurve::new(edge.curve, policy, eligibility.edge_eligible[e]),
            })
            .collect(),
        faces: shell
            .faces
            .into_iter()
            .enumerate()
            .map(|(f, face)| {
                let closure = face
                    .provenance
                    .definition_id
                    .and_then(|id| closure_map.get(&id.get()).copied())
                    .or_else(|| {
                        face.provenance
                            .surface_id
                            .and_then(|id| closure_map.get(&id.get()).copied())
                    });
                CompressedFace {
                    boundaries: face.boundaries,
                    orientation: face.orientation,
                    surface: PolicySurface::with_closure(
                        face.surface,
                        policy,
                        eligibility.face_interior_eligible[f],
                        closure,
                    ),
                    provenance: face.provenance,
                }
            })
            .collect(),
        // The geometric uncertainty is a property of the source representation,
        // not of the wrapped geometry; it must survive the wrap so the
        // tessellator can judge source incidence under the source's own
        // tolerance.
        source_geometric_uncertainty: shell.source_geometric_uncertainty,
    }
}

#[cfg(test)]
mod quotient_tests {
    use super::*;
    use truck_geometry::prelude::{BSplineSurface, KnotVec, Point3};
    use truck_meshalgo::prelude::{InnerSpace, ParametricSurface};

    #[test]
    fn torus_and_sphere_boundaries_keep_038_sampling() {
        use truck_geometry::prelude::{Matrix4, Sphere, Torus};
        use truck_stepio::r#in::step_geometry::{
            Processor, Sphere as StepSphere, TrimmedCurve, UnitCircle,
        };
        let circle = Curve3D::Conic(Conic3D::Circle(Processor::with_transform(
            TrimmedCurve::new(UnitCircle::new(), (0.0, std::f64::consts::TAU)),
            Matrix4::from_scale(4.25),
        )));
        let torus = Surface::ElementarySurface(ElementarySurface::ToroidalSurface(Processor::new(
            Torus::new(Point3::new(0.0, 0.0, 0.0), 3.75, 0.5),
        )));
        let mut shell = CompressedShell {
            vertices: vec![Point3::new(4.25, 0.0, 0.0)],
            edges: vec![CompressedEdge {
                vertices: (0, 0),
                curve: circle,
            }],
            faces: vec![CompressedFace {
                boundaries: vec![vec![(0, true).into()]],
                orientation: true,
                surface: torus,
                provenance: Default::default(),
            }],
            source_geometric_uncertainty: None,
        };
        let eligibility = Eligibility::compute(&shell);
        assert!(!eligibility.edge_eligible[0]);
        assert!(!eligibility.face_interior_eligible[0]);
        let wrapped = wrap_shell(shell.clone(), MeshingPolicy::DEFAULT);
        let range = (0.0, std::f64::consts::TAU);
        assert_eq!(
            wrapped.edges[0].curve.parameter_division(range, 1.0),
            shell.edges[0].curve.parameter_division(range, 1.0),
        );
        let range = (range, range);
        assert_eq!(
            wrapped.faces[0].surface.parameter_division(range, 1.0),
            shell.faces[0].surface.parameter_division(range, 1.0),
        );
        // Source provenance enables the floor for a newly admitted analytic
        // shell, preserving its canonical shared edges.
        shell.faces[0].provenance.surface_id =
            Some(truck_topology::compress::SourceEntityId::new(7));
        let signed = wrap_shell_with_source_metadata(
            shell.clone(),
            MeshingPolicy::DEFAULT,
            &HashMap::new(),
            &HashSet::from([7]),
        );
        let (parameters, _) = signed.edges[0].curve.parameter_division(range.0, 1.0);
        assert!(parameters.len() - 1 >= MeshingPolicy::DEFAULT.angular_floor_segments());
        let (u, v) = signed.faces[0].surface.parameter_division(range, 1.0);
        assert!(u.len() >= MeshingPolicy::DEFAULT.angular_floor_segments());
        assert!(v.len() >= MeshingPolicy::DEFAULT.angular_floor_segments());
        // A positive companion torus in this recovered shell uses the same
        // floor. Without a source-signed carrier the identical shell retains
        // baseline sampling.
        let mut companion = shell.faces[0].clone();
        companion.provenance.surface_id = Some(truck_topology::compress::SourceEntityId::new(8));
        shell.faces.push(companion);
        assert!(
            Eligibility::with_signed_tori(&shell, &HashSet::from([7]))
                .face_interior_eligible
                .iter()
                .all(|eligible| *eligible)
        );
        assert!(
            Eligibility::compute(&shell)
                .face_interior_eligible
                .iter()
                .all(|eligible| !eligible)
        );
        shell.faces[0].surface = Surface::ElementarySurface(ElementarySurface::Sphere(
            Processor::new(StepSphere(Sphere::new(Point3::new(0.0, 0.0, 0.0), 4.25))),
        ));
        assert!(!Eligibility::compute(&shell).edge_eligible[0]);
        let mixed = Eligibility::with_signed_tori(&shell, &HashSet::from([8]));
        assert!(!mixed.edge_eligible[0]);
        assert!(!mixed.face_interior_eligible[1]);
        shell.faces[0].surface = closed_v_spline();
        let mixed_spline = wrap_shell_with_source_metadata(
            shell.clone(),
            MeshingPolicy::DEFAULT,
            &HashMap::new(),
            &HashSet::from([8]),
        );
        assert_eq!(
            mixed_spline.edges[0].curve.parameter_division(range.0, 1.0),
            shell.edges[0].curve.parameter_division(range.0, 1.0)
        );
        assert_eq!(
            mixed_spline.faces[1].surface.parameter_division(range, 1.0),
            shell.faces[1].surface.parameter_division(range, 1.0)
        );
    }

    /// A degree-2 clamped spline over `[0,1]^2` whose `v` direction is a
    /// genuinely closed strip (last control row == first, penultimate row the
    /// reflection `2*P0 - P1`), so a source-closed `v` declaration certifies a
    /// native interval `[0,1]` with period `1.0`.
    fn closed_v_spline() -> Surface {
        let knots = (KnotVec::uniform_knot(2, 2), KnotVec::uniform_knot(2, 2));
        let column = |r: f64, z: f64| {
            let p0 = Point3::new(r, 0.0, z);
            let p1 = Point3::new(r * 0.3, r, z);
            let p2 = Point3::new(p0.x * 2.0 - p1.x, p0.y * 2.0 - p1.y, z);
            vec![p0, p1, p2, p0]
        };
        let control_points = vec![
            column(1.0, 0.0),
            column(1.3, 0.5),
            column(1.0, 1.0),
            column(0.7, 1.5),
        ];
        Surface::BSplineSurface(BSplineSurface::new(knots, control_points))
    }

    #[test]
    fn ordinary_spline_inverse_prefers_native_root_over_extrapolated_root() {
        // x(u) = u - 0.1*u^3 is regular on [0,1], but its extrapolation
        // returns to the same point beyond u=3. The hint must not turn that
        // exterior root into a chart spanning the extrapolated carrier.
        let knots = (KnotVec::bezier_knot(3), KnotVec::bezier_knot(1));
        let controls = [0.0, 1.0 / 3.0, 2.0 / 3.0, 0.9]
            .into_iter()
            .map(|x| vec![Point3::new(x, 0.0, 0.0), Point3::new(x, 1.0, 0.0)])
            .collect();
        let inner = Surface::BSplineSurface(BSplineSurface::new(knots, controls));
        let point = inner.subs(0.25, 0.5);
        let wrapped = PolicySurface::with_closure(
            inner,
            MeshingPolicy::DEFAULT,
            false,
            Some(SplineAxisClosure::OPEN),
        );
        let uv = wrapped
            .search_parameter(point, (3.0, 0.5), 100)
            .expect("an in-domain root exists");
        assert!((0.0..=1.0).contains(&uv.0), "escaped native domain: {uv:?}");
        assert!((wrapped.subs(uv.0, uv.1) - point).magnitude() < 1.0e-6);
    }

    #[test]
    fn ordinary_spline_inverse_finds_native_cell_on_multi_turn_strip() {
        // Authored eight-turn helical strip: a translation by 0.1*u in z
        // gives each physical point an exact exterior representative one
        // turn earlier. The clamped old hint and the domain midpoint both
        // find another turn; the native inverse lies in a remote knot cell.
        let k = 4.0 * (2.0_f64.sqrt() - 1.0) / 3.0;
        let quarters = [
            [(1.0, 0.0), (1.0, k), (k, 1.0), (0.0, 1.0)],
            [(0.0, 1.0), (-k, 1.0), (-1.0, k), (-1.0, 0.0)],
            [(-1.0, 0.0), (-1.0, -k), (-k, -1.0), (0.0, -1.0)],
            [(0.0, -1.0), (k, -1.0), (1.0, -k), (1.0, 0.0)],
        ];
        let mut points = Vec::new();
        for segment in 0..32 {
            for j in 0..3 {
                let (x, y) = quarters[segment % 4][j];
                points.push(Point3::new(x, y, (segment as f64 + j as f64 / 3.0) / 4.0));
            }
        }
        points.push(Point3::new(1.0, 0.0, 8.0));
        let mut knots = vec![0.0; 4];
        for segment in 1..32 {
            knots.extend([segment as f64 / 32.0; 3]);
        }
        knots.extend([1.0; 4]);
        let rows = [0.0, 0.1]
            .into_iter()
            .map(|offset| {
                points
                    .iter()
                    .map(|p| Point3::new(p.x, p.y, p.z + offset))
                    .collect()
            })
            .collect();
        let wrapped = PolicySurface::with_closure(
            Surface::BSplineSurface(BSplineSurface::new(
                (KnotVec::bezier_knot(1), KnotVec::from(knots)),
                rows,
            )),
            MeshingPolicy::DEFAULT,
            false,
            Some(SplineAxisClosure::OPEN),
        );
        let native = (0.4, 0.82);
        let exterior = (10.4, native.1 - 0.125);
        let point = wrapped.subs(native.0, native.1);
        assert!((wrapped.subs(exterior.0, exterior.1) - point).magnitude() < 1.0e-10);
        assert_eq!(
            wrapped.prefer_native_inverse(point, exterior, 100),
            exterior,
            "a face without a topology refusal changed its remote chart"
        );
        let wrapped = wrapped
            .native_inverse_retry()
            .expect("ordinary spline retry");
        let actual = wrapped.prefer_native_inverse(point, exterior, 100);
        assert!(
            wrapped.inverse_is_in_native_domain(actual),
            "native inverse missed: {actual:?}"
        );
        assert!((wrapped.subs(actual.0, actual.1) - point).magnitude() < 1.0e-6);
    }

    #[test]
    fn ordinary_spline_remote_native_root_preserves_sheet_orientation() {
        // The same physical point lies on two sheets of an authored folded
        // helical strip. A remote native root reverses the regular carrier
        // normal. A broad inverse fallback must preserve the original sheet
        // when it lacks topological authority to select that other sheet.
        let k = 4.0 * (2.0_f64.sqrt() - 1.0) / 3.0;
        let quarters = [
            [(1.0, 0.0), (1.0, k), (k, 1.0), (0.0, 1.0)],
            [(0.0, 1.0), (-k, 1.0), (-1.0, k), (-1.0, 0.0)],
            [(-1.0, 0.0), (-1.0, -k), (-k, -1.0), (0.0, -1.0)],
            [(0.0, -1.0), (k, -1.0), (1.0, -k), (1.0, 0.0)],
        ];
        let mut points = Vec::new();
        for segment in 0..32 {
            for j in 0..3 {
                let (x, y) = quarters[segment % 4][j];
                points.push(Point3::new(x, y, (segment as f64 + j as f64 / 3.0) / 4.0));
            }
        }
        points.push(Point3::new(1.0, 0.0, 8.0));
        let mut knots = vec![0.0; 4];
        for segment in 1..32 {
            knots.extend([segment as f64 / 32.0; 3]);
        }
        knots.extend([1.0; 4]);
        let rows = [0.0, 0.0, 1.0]
            .into_iter()
            .map(|offset| {
                points
                    .iter()
                    .map(|p| Point3::new(p.x, p.y, p.z + offset))
                    .collect()
            })
            .collect();
        let wrapped = PolicySurface::with_closure(
            Surface::BSplineSurface(BSplineSurface::new(
                (KnotVec::bezier_knot(2), KnotVec::from(knots)),
                rows,
            )),
            MeshingPolicy::DEFAULT,
            false,
            Some(SplineAxisClosure::OPEN),
        );
        let native = (0.4, 0.82);
        let exterior = (-(1.0_f64 + native.0 * native.0).sqrt(), native.1 - 0.125);
        let point = wrapped.subs(native.0, native.1);
        assert!((wrapped.subs(exterior.0, exterior.1) - point).magnitude() < 1.0e-10);
        let wrapped = wrapped
            .native_inverse_retry()
            .expect("ordinary spline retry");
        let actual = wrapped.prefer_native_inverse(point, exterior, 100);
        assert_eq!(
            actual, exterior,
            "a remote root changed the regular sheet orientation"
        );
    }

    #[test]
    fn ordinary_spline_preserves_exterior_inverse_when_no_native_root_exists() {
        let knots = (KnotVec::bezier_knot(3), KnotVec::bezier_knot(1));
        let controls = [0.0, 1.0 / 3.0, 2.0 / 3.0, 0.9]
            .into_iter()
            .map(|x| vec![Point3::new(x, 0.0, 0.0), Point3::new(x, 1.0, 0.0)])
            .collect();
        let inner = Surface::BSplineSurface(BSplineSurface::new(knots, controls));
        // This point lies beyond the native patch's x<=0.9 control hull.
        let point = inner.subs(2.0, 0.5);
        let wrapped = PolicySurface::with_closure(
            inner,
            MeshingPolicy::DEFAULT,
            false,
            Some(SplineAxisClosure::OPEN),
        );
        let uv = wrapped
            .search_parameter(point, (2.0, 0.5), 100)
            .expect("retain the source's exterior inverse when no native root exists");
        assert!(uv.0 > 1.0);
        assert!((wrapped.subs(uv.0, uv.1) - point).magnitude() < 1.0e-6);
    }

    #[test]
    fn ordinary_spline_preserves_nearly_collapsed_edge_inverse_representative() {
        // S(u,v)=(u,u*v,1e-10*(v-0.1*v^3)). At u=0 the edge is
        // numerically collapsed; distinct v roots cannot select a trim branch
        // reliably. Its exterior inverse must keep the source representative.
        let controls = [0.0, 1.0 / 3.0, 2.0 / 3.0, 0.9];
        for scale in [1.0, 1000.0] {
            let rows = [0.0, 1.0]
                .into_iter()
                .map(|u| {
                    controls
                        .iter()
                        .enumerate()
                        .map(|(j, z)| {
                            Point3::new(
                                scale * u,
                                scale * u * (j as f64 / 3.0),
                                scale * 1.0e-10 * z,
                            )
                        })
                        .collect()
                })
                .collect();
            let inner = Surface::BSplineSurface(BSplineSurface::new(
                (KnotVec::bezier_knot(1), KnotVec::bezier_knot(3)),
                rows,
            ));
            let wrapped = PolicySurface::with_closure(
                inner,
                MeshingPolicy::DEFAULT,
                false,
                Some(SplineAxisClosure::OPEN),
            );
            let source = (0.0, 3.0);
            let point = wrapped.subs(source.0, source.1);
            let actual = wrapped.prefer_native_inverse(point, source, 100);
            assert_eq!(
                actual, source,
                "an ambiguous near-pole root must preserve trim continuity"
            );
        }
    }

    fn v_closed() -> SplineAxisClosure {
        SplineAxisClosure {
            u_closed: false,
            v_closed: true,
        }
    }

    /// The quotient map arithmetic on a native interval that does not start at
    /// zero: `[2.0, 5.0]`, `P = 3.0`.
    #[test]
    fn quotient_map_with_native_interval_not_starting_at_zero() {
        let q = QuotientAxis {
            interval: (2.0, 5.0),
            period: 3.0,
        };
        let cases = [
            (2.0, 2.0),  // native start
            (5.0, 2.0),  // native end == native start (seam)
            (1.5, 4.5),  // one period below
            (5.5, 2.5),  // just above the seam
            (8.5, 2.5),  // two periods above
            (-1.0, 2.0), // exactly a whole period below native start
            (0.0, 3.0),  // mid-deck below
        ];
        for (cover, expected) in cases {
            let got = q.map(cover);
            assert!(
                (got - expected).abs() < 1e-12,
                "map({cover}) = {got}, expected {expected}"
            );
        }
    }

    /// A non-periodic axis passes through exactly: no quotient, byte-identical
    /// evaluation.
    #[test]
    fn non_periodic_axis_passes_through_exactly() {
        let inner = closed_v_spline();
        let wrapped = PolicySurface::with_closure(
            inner.clone(),
            MeshingPolicy::DEFAULT,
            false,
            Some(SplineAxisClosure::OPEN),
        );
        assert_eq!(wrapped.u_quotient, None);
        assert_eq!(wrapped.v_quotient, None);
        for (u, v) in [(0.1, 0.2), (0.5, 0.9), (1.0, 0.0)] {
            let a = wrapped.subs(u, v);
            let b = inner.subs(u, v);
            assert_eq!(a, b, "non-periodic pass-through must be exact");
        }
    }

    /// The certified V axis maps cover UV to the native representative, while
    /// the topology-facing accessors stay unwrapped.
    #[test]
    fn certified_v_axis_maps_evaluation_but_keeps_topology_unwrapped() {
        let inner = closed_v_spline();
        let wrapped = PolicySurface::with_closure(
            inner.clone(),
            MeshingPolicy::DEFAULT,
            false,
            Some(v_closed()),
        );
        let q = wrapped.v_quotient.expect("a certified V quotient");
        assert_eq!(q.interval, (0.0, 1.0));
        assert_eq!(q.period, 1.0);

        // Topology-facing accessors are untouched.
        assert_eq!(wrapped.try_range_tuple(), inner.try_range_tuple());
        assert_eq!(wrapped.u_period(), inner.u_period());
        assert_eq!(wrapped.v_period(), inner.v_period());

        // Evaluator forwarding maps cover -> native on the certified axis.
        let cover = [
            (0.25, 1.7),
            (0.25, -0.25),
            (0.25, 2.25),
            (0.25, 0.7),
            (0.25, 0.0),
        ];
        for (u, v) in cover {
            let expected_v = q.map(v);
            let a = wrapped.subs(u, v);
            let b = inner.subs(u, expected_v);
            assert!(
                (a - b).magnitude() < 1e-12,
                "subs({u},{v}) mapped to v={expected_v} differs from native evaluation"
            );
            let d = wrapped.vder(u, v) - inner.vder(u, expected_v);
            assert!(
                d.magnitude() < 1e-9,
                "vder({u},{v}) mapped to v={expected_v} differs from native evaluation"
            );
        }

        // The analytic-family normal forwarding maps too.
        let n_wrapped = wrapped.normal(0.25, 1.7);
        let n_native = inner.normal(0.25, 0.7);
        assert!(
            (n_wrapped - n_native).magnitude() < 1e-9,
            "normal must evaluate on the native representative"
        );
    }

    /// Stage E (narrowed): an *ordinary* B-spline with no certified periodic
    /// quotient preserves legacy inverse semantics — the strict range rule is a
    /// periodic-policy instrument and must not reject its pre-existing inverse
    /// results (the NIST #1167-wave regression witness). A certified-closed V
    /// axis allows deck-equivalent values while the ordinary U companion stays
    /// bounded: the canonical failure is a non-periodic U search returning
    /// `u = 3.63` for a point whose true `u = 0.765` on a `[0,1]` native domain.
    #[test]
    fn non_periodic_spline_preserves_legacy_inverse_and_certified_axis_stays_bounded() {
        let inner = closed_v_spline();
        // Source open on both axes: no certified period anywhere.
        let open = PolicySurface::with_closure(
            inner.clone(),
            MeshingPolicy::DEFAULT,
            false,
            Some(SplineAxisClosure::OPEN),
        );
        assert!(
            open.accept_inverse_result((3.63, 0.5)),
            "an ordinary spline with no certified quotient must preserve legacy \
             inverse semantics (the exterior root is accepted, as before A-F)"
        );
        assert!(open.accept_inverse_result((0.765, 0.5)));
        // Endpoint values within the justified numerical tolerance are accepted.
        assert!(open.accept_inverse_result((0.0, 0.5)));
        assert!(open.accept_inverse_result((1.0, 0.5)));
        assert!(open.accept_inverse_result((-1.0e-7, 0.5)));
        assert!(open.accept_inverse_result((1.0 + 1.0e-7, 0.5)));

        // The certified V axis allows deck-equivalent cover values, but the
        // non-certified U axis is still bounded.
        let certified =
            PolicySurface::with_closure(inner, MeshingPolicy::DEFAULT, false, Some(v_closed()));
        assert!(
            certified.accept_inverse_result((0.5, 3.63)),
            "a certified periodic axis must allow deck-equivalent values"
        );
        assert!(
            !certified.accept_inverse_result((3.63, 0.5)),
            "the non-certified U axis must still reject an exterior root on a \
             source-certified periodic surface"
        );
    }

    /// Stage E behavioral check: the wrapper's `search_parameter` forwards
    /// through the acceptance gate. On an ordinary spline (no certified
    /// quotient) a legal interior root passes through untouched; on a
    /// source-certified periodic surface the ordinary companion U axis still
    /// rejects an exterior escape while the certified axis admits deck
    /// representatives.
    #[test]
    fn search_parameter_forwards_through_the_acceptance_gate() {
        use truck_meshalgo::prelude::SearchParameter;
        let inner = closed_v_spline();
        let open = PolicySurface::with_closure(
            inner.clone(),
            MeshingPolicy::DEFAULT,
            false,
            Some(SplineAxisClosure::OPEN),
        );
        // A point on the surface at a legal interior UV must be recoverable in
        // domain (the gate must not reject legitimate roots).
        let legal = inner.subs(0.765, 0.5);
        let found = open.search_parameter(legal, None::<(f64, f64)>, 500);
        assert!(found.is_some(), "a legal interior root must be recovered");
        let found = found.unwrap();
        assert!(open.accept_inverse_result(found));
        // The recovered root must be a legal in-domain representative that
        // realizes the point (the B-spline inverse is not guaranteed to return
        // the exact seed, but it must not escape the native domain).
        let realized = inner.subs(found.0, found.1);
        assert!(
            (realized - legal).magnitude() < 1e-2,
            "the accepted root must realize the projected point, got {found:?}"
        );

        // The certified surface keeps the same recovery on the legal root, and
        // the wrapper's inverse must not admit an exterior U escape for a
        // point whose true representative is in-domain.
        let certified =
            PolicySurface::with_closure(inner, MeshingPolicy::DEFAULT, false, Some(v_closed()));
        let certified_found = certified.search_parameter(legal, None::<(f64, f64)>, 500);
        assert!(
            certified_found.is_some(),
            "a legal interior root must be recovered on the certified surface"
        );
        let certified_found = certified_found.unwrap();
        assert!(
            certified.accept_inverse_result(certified_found),
            "the certified-surface inverse result must pass the acceptance gate"
        );
        assert!(
            !certified.accept_inverse_result((3.63, 0.5)),
            "the ordinary companion U axis of a certified periodic surface must \
             reject an exterior root"
        );
        assert!(
            certified.accept_inverse_result((0.765, 3.63)),
            "the certified V axis must admit deck-equivalent cover values"
        );
    }

    /// Stage F theorem: tessellation complexity is invariant to whole-period
    /// deck translation. Equivalent cover intervals of the same physical patch
    /// — `[0.5, 1.5]`, `[1.5, 2.5]` (+P), `[-0.5, 0.5]` (-P), `[0.0, 1.0]`
    /// (pure fundamental domain) — must produce essentially the same
    /// subdivision on the certified axis.
    #[test]
    fn certified_axis_subdivision_is_deck_translation_invariant() {
        use truck_meshalgo::prelude::ParameterDivision2D;
        let inner = closed_v_spline();
        let wrapped =
            PolicySurface::with_closure(inner, MeshingPolicy::DEFAULT, false, Some(v_closed()));
        let tol = 0.05;
        let ranges = [
            ((0.0, 1.0), (0.5, 1.5)),
            ((0.0, 1.0), (1.5, 2.5)),
            ((0.0, 1.0), (-0.5, 0.5)),
            ((0.0, 1.0), (2.5, 3.5)),
            ((0.0, 1.0), (0.0, 1.0)),
        ];
        let counts: Vec<usize> = ranges
            .iter()
            .map(|&range| wrapped.parameter_division(range, tol).1.len())
            .collect();
        let reference = counts[0];
        for (i, &count) in counts.iter().enumerate() {
            assert!(
                count == reference || (count as i64 - reference as i64).unsigned_abs() <= 2,
                "range {i} of {ranges:?} produced {count} V divisions, expected ~{reference} (deck translation must not change subdivision density)"
            );
        }
    }

    /// The same physical patch must also produce the same *physical* grid: the
    /// set of native representatives of a translated cover division is the same.
    #[test]
    fn certified_axis_subdivision_maps_back_to_the_same_physical_grid() {
        use truck_meshalgo::prelude::ParameterDivision2D;
        let inner = closed_v_spline();
        let wrapped =
            PolicySurface::with_closure(inner, MeshingPolicy::DEFAULT, false, Some(v_closed()));
        let q = wrapped.v_quotient.unwrap();
        let tol = 0.05;
        let vdiv_a = wrapped.parameter_division(((0.0, 1.0), (0.5, 1.5)), tol).1;
        let vdiv_b = wrapped.parameter_division(((0.0, 1.0), (1.5, 2.5)), tol).1;
        // Each cover division maps to the same set of native representatives.
        let native_of = |div: &[f64]| {
            let mut set: Vec<f64> = div.iter().map(|&x| q.map(x)).collect();
            set.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            set.dedup_by(|a, b| (*a - *b).abs() < 1.0e-9);
            set
        };
        assert_eq!(native_of(&vdiv_a), native_of(&vdiv_b));
    }

    /// §18 control: analytic periodic surfaces (cylinder, cone, sphere, torus)
    /// never gain a quotient from a source closure — their evaluators are
    /// globally periodic and the Stage C/E/F semantics must not activate on
    /// them. `with_closure` on a non-spline must be a no-op for the quotient.
    #[test]
    fn analytic_periodic_surfaces_never_gain_a_quotient() {
        use truck_geometry::prelude::{Line, Sphere, Torus, Vector3};
        use truck_meshalgo::prelude::{EuclideanSpace, Point3};
        use truck_stepio::r#in::step_geometry::{Processor, Sphere as StepSphere};
        let v_closed = Some(v_closed());
        let cylinder = Processor::new(truck_geometry::prelude::RevolutedCurve::by_revolution(
            Line(Point3::new(10.0, 0.0, 0.0), Point3::new(10.0, 0.0, 10.0)),
            Point3::origin(),
            Vector3::unit_z(),
        ));
        let cone = Processor::new(truck_geometry::prelude::RevolutedCurve::by_revolution(
            Line(Point3::new(10.0, 0.0, 0.0), Point3::new(0.0, 0.0, 10.0)),
            Point3::origin(),
            Vector3::unit_z(),
        ));
        let surfaces = [
            Surface::ElementarySurface(ElementarySurface::CylindricalSurface(cylinder)),
            Surface::ElementarySurface(ElementarySurface::ConicalSurface(cone)),
            Surface::ElementarySurface(ElementarySurface::Sphere(Processor::new(StepSphere(
                Sphere::new(Point3::origin(), 5.0),
            )))),
            Surface::ElementarySurface(ElementarySurface::ToroidalSurface(Processor::new(
                Torus::new(Point3::origin(), 8.0, 2.0),
            ))),
        ];
        for surface in &surfaces {
            let wrapped = PolicySurface::with_closure(
                surface.clone(),
                MeshingPolicy::DEFAULT,
                false,
                v_closed,
            );
            assert_eq!(
                wrapped.u_quotient, None,
                "an analytic surface must never gain a U quotient"
            );
            assert_eq!(
                wrapped.v_quotient, None,
                "an analytic surface must never gain a V quotient"
            );
            // Evaluation forwards byte-identically to the unwrapped path.
            for (u, v) in [(0.3, 1.0), (1.7, -0.5), (2.0, 4.0)] {
                assert_eq!(wrapped.subs(u, v), surface.subs(u, v));
            }
            // The inverse gate does not bound analytic axes either.
            assert!(
                wrapped.accept_inverse_result((3.63, 0.5)),
                "analytic axes are not bounded by the inverse gate"
            );
        }
    }
}
