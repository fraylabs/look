use proptest::{prelude::*, property_test};
use truck_geometry::prelude::*;

#[property_test]
fn sp_test(#[strategy = -10f64..=10f64] t: f64) {
    let curve = UnitHyperbola::<Point2>::new();
    let p = curve.subs(t);
    prop_assert_near!(curve.search_parameter(p, None, 0).unwrap(), t);
}

#[test]
fn snp_test() {
    let (t, r) = (2.0, 1.0);
    let curve = UnitHyperbola::<Point2>::new();
    let p = curve.subs(t);
    let q = p + r * Vector2::new(-p.x, p.y);
    let t = curve.search_nearest_parameter(q, None, 0).unwrap();
    let p = curve.subs(t);
    let dot = curve.der(t).dot(q - p);
    assert!(dot.so_small(), "{t} {dot}");
}

#[test]
fn sp_negative_test() {
    let curve = UnitHyperbola::<Point2>::new();
    let q = Point2::new(-1.0, 0.0);
    assert!(curve.search_parameter(q, None, 0).is_none());
}

// Normalized endpoints from FreeCAD hyperbola #6174 in the real STEP
// face-loss repro. Quartic roots are sinh(t), not curve parameters.
#[test]
fn nearest_parameter_on_negative_branch_step_endpoints() {
    let curve = UnitHyperbola::<Point2>::new();
    for t in [-2.761967, -2.742247, 2.761967, 2.742247] {
        let point = curve.subs(t);
        let nearest = curve.search_nearest_parameter(point, None, 0).unwrap();
        assert!((nearest - t).abs() < 1.0e-8, "expected {t}, got {nearest}");
        assert!(point.distance(curve.subs(nearest)) < 1.0e-7);
    }
}

#[test]
fn nearest_parameter_refuses_nonfinite_points() {
    let curve = UnitHyperbola::<Point2>::new();
    for point in [Point2::new(f64::NAN, 0.0), Point2::new(1.0, f64::INFINITY)] {
        assert!(curve.search_nearest_parameter(point, None, 0).is_none());
    }
}
