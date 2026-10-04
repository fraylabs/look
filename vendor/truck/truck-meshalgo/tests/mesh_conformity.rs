//! Shared BREP edges must produce the same mesh boundary on every incident face.
use truck_meshalgo::prelude::*;
use truck_modeling::*;
use truck_topology::shell::ShellCondition;

fn assert_closed(solid: &Solid) {
    let shell = solid.boundaries()[0].compress();
    let meshed = shell.robust_triangulation(0.01);
    assert!(
        meshed.faces.iter().all(|face| face
            .surface
            .as_ref()
            .is_some_and(|mesh| !mesh.faces().is_empty())),
        "every BREP face must survive"
    );
    let mut mesh = meshed.to_polygon();
    mesh.put_together_same_attrs(TOLERANCE)
        .remove_unused_attrs();
    assert_eq!(
        mesh.shell_condition(),
        ShellCondition::Closed,
        "every mesh edge must have exactly two opposite incident triangles"
    );
}

#[test]
fn cylinder_seam_is_closed_and_manifold() {
    let vertex = builder::vertex(Point3::new(1.0, 0.0, 0.0));
    let circle: Wire = builder::rsweep(
        &vertex,
        Point3::origin(),
        Vector3::unit_z(),
        Rad(std::f64::consts::TAU),
        2,
    );
    let disk: Face = builder::try_attach_plane(vec![circle]).expect("planar circle"); // H-1: test fixture
    let cylinder: Solid = builder::tsweep(&disk, Vector3::new(0.0, 0.0, 2.0));
    assert_closed(&cylinder);
}

fn filleted_box() -> Solid {
    let points = [
        (-1.0, -2.0),
        (1.0, -2.0),
        (2.0, -1.0),
        (2.0, 1.0),
        (1.0, 2.0),
        (-1.0, 2.0),
        (-2.0, 1.0),
        (-2.0, -1.0),
    ];
    let vertices = builder::vertices(points.map(|(x, y)| Point3::new(x, y, 0.0)));
    let r = 1.0 / 2.0_f64.sqrt();
    let arcs = [
        (1.0 + r, -1.0 - r),
        (1.0 + r, 1.0 + r),
        (-1.0 - r, 1.0 + r),
        (-1.0 - r, -1.0 - r),
    ];
    let mut edges: Vec<Edge> = Vec::new();
    for i in 0..8 {
        let a = &vertices[i];
        let b = &vertices[(i + 1) % 8];
        if i % 2 == 0 {
            edges.push(builder::line(a, b));
        } else {
            let (x, y) = arcs[i / 2];
            edges.push(builder::circle_arc(a, b, Point3::new(x, y, 0.0)));
        }
    }
    let disk: Face =
        builder::try_attach_plane(vec![Wire::from(edges)]).expect("planar rounded rectangle"); // H-1: test fixture
    builder::tsweep(&disk, Vector3::new(0.0, 0.0, 3.0))
}

#[test]
fn filleted_box_is_closed_and_manifold() {
    assert_closed(&filleted_box());
}

#[test]
fn shell_endpoints_share_source_vertex() {
    let solid = filleted_box();
    let mut shell = solid.boundaries()[0].compress();
    for (i, edge) in shell.edges.iter_mut().enumerate() {
        if let Curve::Line(line) = &mut edge.curve {
            let delta = Vector3::new((i + 1) as f64 * TOLERANCE * 4.0, 0.0, 0.0);
            line.0 += delta;
            line.1 += delta;
        }
    }
    let meshed = shell.robust_triangulation(0.01);
    assert!(
        meshed.faces.iter().all(|face| face.surface.is_some()),
        "tolerance-compatible endpoints must mesh"
    );
    for edge in &meshed.edges {
        assert_eq!(
            edge.curve.first().copied(),
            shell.vertices.get(edge.vertices.0),
            "edge starts at its shared BREP vertex"
        );
        assert_eq!(
            edge.curve.last().copied(),
            shell.vertices.get(edge.vertices.1),
            "edge ends at its shared BREP vertex"
        );
    }
    let mut mesh = meshed.to_polygon();
    mesh.put_together_same_attrs(TOLERANCE)
        .remove_unused_attrs();
    assert_eq!(mesh.shell_condition(), ShellCondition::Closed);
}
