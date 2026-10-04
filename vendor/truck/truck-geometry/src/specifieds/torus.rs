use super::*;
use std::f64::consts::PI;

impl Torus {
    /// constructor
    #[inline(always)]
    pub fn new(center: Point3, large_radius: f64, small_radius: f64) -> Self {
        if large_radius <= 0.0 || small_radius <= 0.0 {
            panic!("radius must be larger than 0");
        }
        Self {
            center,
            large_radius,
            small_radius,
        }
    }

    /// The two meridian sections at a point's azimuth. A spindle torus's
    /// folded sheet has negative signed radial coordinate, so its parameter u
    /// differs from the point's azimuth by pi. Both candidates are needed even
    /// when the observed radial distance exceeds the major radius.
    fn meridian_parameters(&self, point: Point3) -> [(f64, f64); 2] {
        let relative = point - self.center;
        let rho = relative.x.hypot(relative.y);
        let azimuth = relative.y.atan2(relative.x);
        [1.0, -1.0].map(|sign| {
            let u = (azimuth + if sign < 0.0 { PI } else { 0.0 }).rem_euclid(2.0 * PI);
            let v = relative
                .z
                .atan2(sign * rho - self.large_radius)
                .rem_euclid(2.0 * PI);
            (u, v)
        })
    }

    /// Sign of the parameterisation's radial Jacobian. The unit meridian
    /// vector points opposite to du cross dv on the folded spindle sheet.
    fn radial_orientation(&self, v: f64) -> f64 {
        let radius = self.large_radius + self.small_radius * v.cos();
        if radius == 0.0 {
            0.0
        } else {
            radius.signum()
        }
    }

    /// get center
    #[inline(always)]
    pub const fn center(&self) -> Point3 {
        self.center
    }

    /// get large radius
    #[inline(always)]
    pub const fn large_radius(&self) -> f64 {
        self.large_radius
    }

    /// get small radius
    #[inline(always)]
    pub const fn small_radius(&self) -> f64 {
        self.small_radius
    }
}

impl ParametricSurface for Torus {
    type Point = Point3;
    type Vector = Vector3;
    #[inline(always)]
    fn der_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Self::Vector {
        let ((su, cu), (sv, cv)) = (u.sin_cos(), v.sin_cos());
        let center = match (m, n) {
            (0, 0) => self.center.to_vec(),
            _ => Vector3::zero(),
        };
        let u_z = if m == 0 { 1.0 } else { 0.0 };
        let u_part = match m % 4 {
            0 => Vector3::new(cu, su, u_z),
            1 => Vector3::new(-su, cu, 0.0),
            2 => Vector3::new(-cu, -su, 0.0),
            _ => Vector3::new(su, -cu, 0.0),
        };
        let r0 = if n == 0 { self.large_radius } else { 0.0 };
        let r1 = self.small_radius;
        let v_part_d2 = match n % 4 {
            0 => Vector2::new(r0 + r1 * cv, r1 * sv),
            1 => Vector2::new(-r1 * sv, r1 * cv),
            2 => Vector2::new(-r1 * cv, -r1 * sv),
            _ => Vector2::new(r1 * sv, -r1 * cv),
        };
        let v_part = Vector3::new(v_part_d2.x, v_part_d2.x, v_part_d2.y);
        center + u_part.mul_element_wise(v_part)
    }
    #[inline(always)]
    fn subs(&self, u: f64, v: f64) -> Point3 {
        let sr = self.small_radius() * Vector2::new(f64::cos(v), f64::sin(v));
        let lr = (self.large_radius() + sr.x) * Vector2::new(f64::cos(u), f64::sin(u));
        self.center() + Vector3::new(lr.x, lr.y, sr.y)
    }
    #[inline(always)]
    fn uder(&self, u: f64, v: f64) -> Vector3 {
        let sr = self.small_radius() * f64::cos(v);
        let lr = (self.large_radius() + sr) * Vector2::new(f64::cos(u), f64::sin(u));
        Vector3::new(-lr.y, lr.x, 0.0)
    }
    #[inline(always)]
    fn vder(&self, u: f64, v: f64) -> Vector3 {
        let sv = self.small_radius() * Vector2::new(-f64::sin(v), f64::cos(v));
        Vector3::new(sv.x * f64::cos(u), sv.x * f64::sin(u), sv.y)
    }
    #[inline(always)]
    fn uuder(&self, u: f64, v: f64) -> Vector3 {
        let sr = self.small_radius() * f64::cos(v);
        let lr = (self.large_radius() + sr) * Vector2::new(f64::cos(u), f64::sin(u));
        Vector3::new(-lr.x, -lr.y, 0.0)
    }
    #[inline(always)]
    fn uvder(&self, u: f64, v: f64) -> Vector3 {
        let sr = -self.small_radius() * f64::sin(v);
        let lr = sr * Vector2::new(f64::cos(u), f64::sin(u));
        Vector3::new(-lr.y, lr.x, 0.0)
    }
    #[inline(always)]
    fn vvder(&self, u: f64, v: f64) -> Vector3 {
        let sv = -self.small_radius() * Vector2::new(f64::cos(v), f64::sin(v));
        Vector3::new(sv.x * f64::cos(u), sv.x * f64::sin(u), sv.y)
    }
    #[inline(always)]
    fn parameter_range(&self) -> (ParameterRange, ParameterRange) {
        const RANGE: (Bound<f64>, Bound<f64>) = (Bound::Included(0.0), Bound::Excluded(2.0 * PI));
        (RANGE, RANGE)
    }
    #[inline(always)]
    fn u_period(&self) -> Option<f64> {
        Some(2.0 * PI)
    }
    #[inline(always)]
    fn v_period(&self) -> Option<f64> {
        Some(2.0 * PI)
    }
}

impl ParametricSurface3D for Torus {
    #[inline(always)]
    fn normal(&self, u: f64, v: f64) -> Vector3 {
        let sv = Vector2::new(f64::cos(v), f64::sin(v));
        Vector3::new(sv.x * f64::cos(u), sv.x * f64::sin(u), sv.y) * self.radial_orientation(v)
    }
    #[inline(always)]
    fn normal_uder(&self, u: f64, v: f64) -> Vector3 {
        let sv = Vector2::new(f64::cos(v), f64::sin(v));
        Vector3::new(-sv.x * f64::sin(u), sv.x * f64::cos(u), 0.0) * self.radial_orientation(v)
    }
    #[inline(always)]
    fn normal_vder(&self, u: f64, v: f64) -> Vector3 {
        let sv = Vector2::new(-f64::sin(v), f64::cos(v));
        Vector3::new(sv.x * f64::cos(u), sv.x * f64::sin(u), sv.y) * self.radial_orientation(v)
    }
}

impl BoundedSurface for Torus {}

impl SearchParameter<D2> for Torus {
    type Point = Point3;
    fn search_parameter<H: Into<SPHint2D>>(
        &self,
        point: Point3,
        _: H,
        _: usize,
    ) -> Option<(f64, f64)> {
        let ctx = ToleranceCtx::unscaled_legacy();
        self.meridian_parameters(point)
            .into_iter()
            .find(|&(u, v)| ctx.near_pt(self.subs(u, v), point))
    }
}

impl SearchNearestParameter<D2> for Torus {
    type Point = Point3;
    fn search_nearest_parameter<H: Into<SPHint2D>>(
        &self,
        point: Point3,
        _: H,
        _: usize,
    ) -> Option<(f64, f64)> {
        let ctx = ToleranceCtx::unscaled_legacy();
        let relative = point - self.center;
        // The axis has no distinguished azimuth. Retain the existing refusal
        // for a nearest-point query at this singular set.
        if ctx.is_small_len(relative.x.hypot(relative.y)) {
            return None;
        }
        self.meridian_parameters(point)
            .into_iter()
            .min_by(|&(u0, v0), &(u1, v1)| {
                self.subs(u0, v0)
                    .distance2(point)
                    .total_cmp(&self.subs(u1, v1).distance2(point))
            })
    }
}

impl ParameterDivision2D for Torus {
    fn parameter_division(
        &self,
        (urange, vrange): ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let circle = UnitCircle::<Point2>::new();
        let utol = tol / (self.small_radius() + self.large_radius());
        let (udiv, _) = circle.parameter_division(urange, utol);
        let vtol = tol / self.small_radius();
        let (vdiv, _) = circle.parameter_division(vrange, vtol);
        (udiv, vdiv)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn torus_normal_derivatives_match_finite_difference_on_both_sheets() {
        let h = 1e-6;
        for torus in [
            Torus::new(Point3::origin(), 3.0, 1.0),
            Torus::new(Point3::origin(), 1.55, 5.0),
        ] {
            for u in [0.0, 0.3, 1.0, 2.0, 5.0] {
                for v in [0.2, 1.0, 2.8, 3.1, 4.0] {
                    let du = (torus.normal(u + h, v) - torus.normal(u - h, v)) / (2.0 * h);
                    let dv = (torus.normal(u, v + h) - torus.normal(u, v - h)) / (2.0 * h);
                    assert!((torus.normal_uder(u, v) - du).magnitude() < 1e-8);
                    assert!((torus.normal_vder(u, v) - dv).magnitude() < 1e-8);
                    assert!(
                        (torus.normal(u, v) - torus.uder(u, v).cross(torus.vder(u, v)).normalize())
                            .magnitude()
                            < 1e-12
                    );
                }
            }
        }
    }
}
