use truck_geometry::prelude::*;
use truck_meshalgo::prelude::*;
use truck_modeling::{Curve, Surface};
use truck_topology::compress::*;

fn planar_shell() -> CompressedShell<Point3, Curve, Surface> {
    let a = Point3::new(1.00005, 0., 0.);
    let b = Point3::new(0., 1.00005, 0.);
    CompressedShell {
        vertices: vec![a, b],
        edges: vec![
            CompressedEdge {
                vertices: (0, 1),
                curve: Curve::Circle(Processor::new(TrimmedCurve::new(
                    UnitCircle::<Point3>::new(),
                    (0., std::f64::consts::FRAC_PI_2),
                ))),
            },
            CompressedEdge {
                vertices: (1, 0),
                curve: Curve::Line(Line(b, a)),
            },
        ],
        faces: vec![CompressedFace {
            boundaries: vec![vec![
                CompressedEdgeIndex {
                    index: 0,
                    orientation: true,
                },
                CompressedEdgeIndex {
                    index: 1,
                    orientation: true,
                },
            ]],
            orientation: true,
            surface: Surface::Plane(Plane::new(
                Point3::origin(),
                Point3::new(1., 0., 0.),
                Point3::new(0., 1., 0.),
            )),
            provenance: Default::default(),
        }],
        source_geometric_uncertainty: None,
    }
}

#[test]
fn approximate_vertices_do_not_deform_the_edge_carrier() {
    let mut shell = planar_shell();
    shell.source_geometric_uncertainty = Some(0.0001);
    let mesh = shell.triangulation(0.001);
    assert!(mesh.faces[0].surface.is_some());
    // Keep the analytic carrier and its interior samples unchanged while
    // realizing admitted endpoints at the shared source vertices. Otherwise
    // this edge and its adjacent line leave cracks in a conforming shell.
    assert_eq!(shell.edges[0].curve.subs(0.), Point3::new(1., 0., 0.));
    let samples = &mesh.edges[0].curve;
    assert_eq!(samples.first(), Some(&shell.vertices[0]));
    assert_eq!(samples.last(), Some(&shell.vertices[1]));
    assert!(samples.len() > 2);
    for point in &samples[1..samples.len() - 1] {
        assert!((point.distance(Point3::origin()) - 1.).abs() < 1e-10);
    }
}

#[test]
fn a_panicking_face_evaluator_does_not_discard_neighboring_faces() {
    let mut shell = planar_shell();
    shell.faces.push(shell.faces[0].clone());
    shell.faces[1].orientation = false;
    let mut poisoned = shell.faces[1].surface.clone();
    poisoned.transform_by(Matrix4::from_translation(Vector3::new(2., 0., 0.)));
    shell.faces[1].surface = poisoned;
    let outcome = shell.robust_triangulation_with_lattice_outcome(0.001, |surface| {
        if surface.subs(0., 0.).x > 1. {
            panic!("unexpected surface evaluator failure");
        }
        unevidenced_lattice(surface)
    });
    assert!(outcome.shell.faces[0].surface.is_some());
    assert!(outcome.shell.faces[1].surface.is_none());
    assert_eq!(
        outcome.face_failures[1].as_ref().unwrap().reason,
        truck_meshalgo::tessellation::TessellationFailureReason::KernelPanicked
    );
}

#[test]
fn endpoints_outside_the_error_budget_still_refuse_the_face() {
    let mut shell = planar_shell();
    shell.vertices[0].x = 1.1;
    shell.vertices[1].y = 1.1;
    shell.edges[1].curve = Curve::Line(Line(shell.vertices[1], shell.vertices[0]));
    let mesh = shell.triangulation(0.001);
    assert!(mesh.faces[0].surface.is_none());
}
