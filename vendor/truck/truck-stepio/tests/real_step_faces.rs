// Source licences for embedded real-model STEP data: real_step_faces.NOTICE.md.
use truck_stepio::r#in::{convert::FaceLossReason, Table};

fn table(data: &str) -> Table {
    use std::str::FromStr;
    let parsed = ruststep::ast::DataSection::from_str(&format!("DATA;{data}ENDSEC;")).unwrap();
    Table::from_owned_data_section(parsed)
}

#[test]
fn freecad_hyperbola_face_survives_conversion() {
    // Geometric reproduction of Apollo FreeCAD face #11360.
    let table = table(
        r#"
#6154=VERTEX_POINT('',#6155);
#6155=CARTESIAN_POINT('',(29.792348169434,15.671104406602,32.986157));
#6162=EDGE_CURVE('',#6163,#6154,#6165,.T.);
#6163=VERTEX_POINT('',#6164);
#6164=CARTESIAN_POINT('',(30.599843085425,15.764582259443,
    33.110600947936));
#6165=ELLIPSE('',#6166,6.42642895366,6.275569927951);
#6166=AXIS2_PLACEMENT_3D('',#6167,#6168,#6169);
#6167=CARTESIAN_POINT('',(29.458515318402,21.938,32.934709644592));
#6168=DIRECTION('',(0.152313001304,0.,-0.98833230729));
#6169=DIRECTION('',(-0.98833230729,-0.,-0.152313001304));
#6171=EDGE_CURVE('',#6172,#6163,#6174,.T.);
#6172=VERTEX_POINT('',#6173);
#6173=CARTESIAN_POINT('',(30.619021281933,15.640620119031,32.986157));
#6174=HYPERBOLA('',#6175,0.817966814837,0.798765221287);
#6175=AXIS2_PLACEMENT_3D('',#6176,#6177,#6178);
#6176=CARTESIAN_POINT('',(29.628920970605,21.938,39.410743980535));
#6177=DIRECTION('',(0.98833230729,-0.,0.152313001304));
#6178=DIRECTION('',(0.152313001304,0.,-0.98833230729));
#11360=ADVANCED_FACE('',(#11361),#11372,.T.);
#11361=FACE_BOUND('',#11362,.T.);
#11362=EDGE_LOOP('',(#11363,#11364,#11371));
#11363=ORIENTED_EDGE('',*,*,#6162,.T.);
#11364=ORIENTED_EDGE('',*,*,#11365,.T.);
#11365=EDGE_CURVE('',#6154,#6172,#11366,.T.);
#11366=CIRCLE('',#11367,6.3);
#11367=AXIS2_PLACEMENT_3D('',#11368,#11369,#11370);
#11368=CARTESIAN_POINT('',(30.437344,21.938,32.986157));
#11369=DIRECTION('',(0.,0.,1.));
#11370=DIRECTION('',(1.,0.,0.));
#11371=ORIENTED_EDGE('',*,*,#6171,.T.);
#11372=CONICAL_SURFACE('',#11373,6.3,0.785398163397);
#11373=AXIS2_PLACEMENT_3D('',#11374,#11375,#11376);
#11374=CARTESIAN_POINT('',(30.437344,21.938,32.986157));
#11375=DIRECTION('',(-0.,-0.,-1.));
#11376=DIRECTION('',(1.,0.,0.));
#25264=OPEN_SHELL('',(#11360));
"#,
    );
    let (shell, losses) = table
        .to_compressed_shell_with_losses(25264, &table.shell[&25264])
        .unwrap();
    assert!(losses.is_empty(), "{losses:?}");
    assert_eq!(shell.faces.len(), 1);
}

#[test]
fn zero_major_torus_is_a_typed_face_refusal() {
    let table = table(
        r#"
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=TOROIDAL_SURFACE('',#4,0.0,5.0);
#6=ADVANCED_FACE('',(),#5,.T.);
#7=OPEN_SHELL('',(#6));
"#,
    );
    let (shell, losses) = table
        .to_compressed_shell_with_losses(7, &table.shell[&7])
        .unwrap();
    assert!(shell.faces.is_empty());
    assert_eq!(losses.len(), 1);
    assert_eq!(losses[0].reason, FaceLossReason::SurfaceConversionFailed);
}

#[test]
fn approximate_edge_incidence_is_checked_at_tessellation_tolerance() {
    use truck_meshalgo::tessellation::MeshableShape;
    // Two circles declared incident to the same vertex, displaced by 5e-5.
    // Conversion must preserve the edges until the mesh tolerance is known.
    let table = table(
        r#"
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=CIRCLE('',#4,1.);
#6=PLANE('',#4);
#7=CARTESIAN_POINT('',(1.00005,0.,0.));
#8=VERTEX_POINT('',#7);
#9=EDGE_CURVE('',#8,#8,#5,.T.);
#10=ORIENTED_EDGE('',*,*,#9,.T.);
#11=EDGE_LOOP('',(#10));
#12=FACE_OUTER_BOUND('',#11,.T.);
#13=ADVANCED_FACE('',(#12),#6,.T.);
#14=OPEN_SHELL('',(#13));
"#,
    );
    let (shell, losses) = table
        .to_compressed_shell_with_losses(14, &table.shell[&14])
        .unwrap();
    assert!(losses.is_empty(), "{losses:?}");
    let mesh = shell.triangulation(0.001);
    assert!(mesh.faces[0].surface.is_some());
}

#[test]
fn reversed_surface_curve_is_inverted_once() {
    use ruststep::tables::EntityTable;
    use truck_stepio::r#in::{step_geometry::*, EdgeCurveHolder};
    let table = table(
        r#"
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=CIRCLE('',#4,1.);
#6=SURFACE_CURVE('',#5,(),.CURVE_3D.);
#7=CARTESIAN_POINT('',(1.,0.,0.));
#8=CARTESIAN_POINT('',(0.,1.,0.));
#9=VERTEX_POINT('',#7);
#10=VERTEX_POINT('',#8);
#11=EDGE_CURVE('',#9,#10,#6,.F.);
"#,
    );
    let edge = EntityTable::<EdgeCurveHolder>::get_owned(&table, 11).unwrap();
    let curve = edge.parse_curve3d().unwrap();
    assert!(curve.front().distance(Point3::new(1., 0., 0.)) < 1e-10);
    assert!(curve.back().distance(Point3::new(0., 1., 0.)) < 1e-10);
    let (a, b) = curve.range_tuple();
    let mid = curve.subs((a + b) / 2.);
    assert!(mid.x < 0. && mid.y < 0., "wrong source arc: {mid:?}");
}

#[test]
fn one_panicking_surface_does_not_discard_neighboring_faces() {
    let table = table(
        r#"
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=BEZIER_SURFACE('',3,3,(),.UNSPECIFIED.,.F.,.F.,.F.);
#6=SPHERICAL_SURFACE('',#4,1.);
#7=ADVANCED_FACE('',(),#5,.T.);
#8=ADVANCED_FACE('',(),#6,.T.);
#9=OPEN_SHELL('',(#7,#8));
"#,
    );
    let (shell, losses) = table
        .to_compressed_shell_with_losses(9, &table.shell[&9])
        .unwrap();
    assert_eq!(shell.faces.len(), 1);
    assert_eq!(losses.len(), 1);
    assert_eq!(losses[0].reason, FaceLossReason::SurfaceConversionFailed);
}

#[test]
fn voron2_small_cylinder_retains_carriers_and_refuses_off_surface_vertices() {
    use truck_geometry::prelude::*;
    use truck_meshalgo::tessellation::MeshableShape;
    let table = table(
        r#"
#101951=CYLINDRICAL_SURFACE('',#1879244,0.005);
#142480=FACE_OUTER_BOUND('',#244816,.T.);
#244816=EDGE_LOOP('',(#1245264,#1245265));
#437651=LINE('',#2774758,#585793);
#437680=LINE('',#2774855,#585822);
#585793=VECTOR('',#2109006,1000.);
#585822=VECTOR('',#2109059,1000.);
#738934=VERTEX_POINT('',#2774755);
#738935=VERTEX_POINT('',#2774757);
#927388=EDGE_CURVE('',#738934,#738935,#437651,.T.);
#927425=EDGE_CURVE('',#738934,#738935,#437680,.T.);
#1245264=ORIENTED_EDGE('',*,*,#927388,.F.);
#1245265=ORIENTED_EDGE('',*,*,#927425,.T.);
#1749736=ADVANCED_FACE('',(#142480),#101951,.T.);
#1879244=AXIS2_PLACEMENT_3D('',#2774854,#2109057,#2109058);
#2109006=DIRECTION('',(-1.,1.11022302462478E-16,-1.22460626911908E-16));
#2109057=DIRECTION('center_axis',(1.,-1.11022302462475E-16,1.22460635382238E-16));
#2109058=DIRECTION('ref_axis',(1.22460635382238E-16,-1.23259516444169E-29,
-1.));
#2109059=DIRECTION('',(-1.,1.11022302462478E-16,-1.22460626911908E-16));
#2774755=CARTESIAN_POINT('',(174.040900232234,-98.5449002322341,-0.0524997677659169));
#2774757=CARTESIAN_POINT('',(173.195899767766,-98.5449002322341,-0.0524997677659163));
#2774758=CARTESIAN_POINT('',(173.1934,-98.5474,-0.0550000000000006));
#2774854=CARTESIAN_POINT('Origin',(173.1934,-98.5424,-0.0549999999999984));
#2774855=CARTESIAN_POINT('',(173.1934,-98.5424,-0.0499999999999989));
#9999999=OPEN_SHELL('',(#1749736));
"#,
    );
    let (mut shell, losses) = table
        .to_compressed_shell_with_losses(9999999, &table.shell[&9999999])
        .unwrap();
    assert!(losses.is_empty(), "{losses:?}");
    assert!(
        shell.edges[0]
            .curve
            .subs(0.5)
            .distance(shell.edges[1].curve.subs(0.5))
            > 0.005
    );
    // The original assembly declares 0.01-unit source uncertainty.
    shell.source_geometric_uncertainty = Some(0.01);
    // Source incidence admits these two distinct line carriers, but their
    // shared vertices are 0.001464 units off the cylinder. Conforming meshes
    // must retain the shared vertices; the source uncertainty does not widen
    // the surface projection's 0.001-unit mesh bound or authorize a new trim.
    for vertex in &shell.vertices {
        let surface = &shell.faces[0].surface;
        let (u, v) = surface.search_nearest_parameter(*vertex, None, 100).unwrap();
        assert!(surface.subs(u, v).distance(*vertex) > 0.001);
    }
    let mesh = shell.triangulation(0.001);
    assert!(mesh.faces[0].surface.is_none());
}

#[test]
fn numerically_incident_line_vertices_remain_exactly_shared() {
    use ruststep::tables::EntityTable;
    use truck_geometry::prelude::*;
    use truck_stepio::r#in::EdgeCurveHolder;
    let table = table(
        r#"
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(1.,0.,0.));
#3=VECTOR('',#2,1.);
#4=LINE('',#1,#3);
#5=CARTESIAN_POINT('',(1.,1.E-13,0.));
#6=CARTESIAN_POINT('',(2.,1.E-13,0.));
#7=VERTEX_POINT('',#5);
#8=VERTEX_POINT('',#6);
#9=EDGE_CURVE('',#7,#8,#4,.T.);
"#,
    );
    let edge = EntityTable::<EdgeCurveHolder>::get_owned(&table, 9).unwrap();
    let curve = edge.parse_curve3d().unwrap();
    assert_eq!(curve.front(), Point3::new(1., 1.0E-13, 0.));
    assert_eq!(curve.back(), Point3::new(2., 1.0E-13, 0.));
}

#[test]
fn negative_major_torus_preserves_occt_geometry_parameterisation() {
    use ruststep::tables::EntityTable;
    use truck_stepio::r#in::{step_geometry::*, ToroidalSurfaceHolder};
    // Synthetic SolidWorks/ProE convention, as implemented by OCCT 7.9.3:
    // StepToGeom takes abs(major), StepToTopoDS reverses the face sense.
    let table = table(
        r#"
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=TOROIDAL_SURFACE('',#4,-1.55,5.0);
"#,
    );
    let source = EntityTable::<ToroidalSurfaceHolder>::get_owned(&table, 5).unwrap();
    let surface = ToroidalSurface::try_from(&source).expect("signed STEP major radius");
    assert_eq!(surface.entity().large_radius(), 1.55);
    assert!(surface.orientation());
    for (u, v) in [(0.4_f64, 0.2_f64), (1.1, 2.7), (2.2, 3.1), (4.1, 3.8)] {
        let radius = 1.55 + 5.0 * v.cos();
        let expected = Point3::new(radius * u.cos(), radius * u.sin(), 5.0 * v.sin());
        // The carrier keeps native STEP UV; face sense is separate.
        assert!(surface.subs(u, v).distance(expected) < 1e-12);
        let uv = surface.search_parameter(expected, None, 100).unwrap();
        assert!(surface.subs(uv.0, uv.1).distance(expected) < 1e-10);
        let nearest = surface
            .search_nearest_parameter(expected, None, 100)
            .unwrap();
        assert!(surface.subs(nearest.0, nearest.1).distance(expected) < 1e-10);
        let outward = Vector3::new(v.cos() * u.cos(), v.cos() * u.sin(), v.sin()) * radius.signum();
        assert!((surface.normal(u, v) - outward).magnitude() < 1e-12);
    }
}

#[test]
fn spindle_torus_folded_sheet_inverse_and_normal_match_its_derivatives() {
    use truck_stepio::r#in::step_geometry::*;
    let torus = Torus::new(Point3::origin(), 1.55, 5.0);
    let (u, v) = (0.7, 3.1);
    let point = torus.subs(u, v);
    let uv = torus
        .search_parameter(point, None, 100)
        .expect("inner spindle sheet inverse");
    assert!(torus.subs(uv.0, uv.1).distance(point) < 1e-10);
    let nearest = torus.search_nearest_parameter(point, None, 100).unwrap();
    assert!(torus.subs(nearest.0, nearest.1).distance(point) < 1e-10);
    assert!(
        (torus.normal(u, v) - torus.uder(u, v).cross(torus.vder(u, v)).normalize()).magnitude()
            < 1e-12
    );
}

#[test]
fn negative_torus_pcurve_keeps_source_uv_coordinates() {
    use ruststep::tables::EntityTable;
    use truck_stepio::r#in::{step_geometry::*, PcurveHolder};
    let table = table(
        r#"
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=TOROIDAL_SURFACE('',#4,-1.55,5.0);
#6=CARTESIAN_POINT('',(0.4,0.2));
#7=DIRECTION('',(1.,0.));
#8=VECTOR('',#7,1.);
#9=LINE('',#6,#8);
#10=REPRESENTATION_CONTEXT('','');
#11=DEFINITIONAL_REPRESENTATION('',(#9),#10);
#12=PCURVE('',#5,#11);
"#,
    );
    let source = EntityTable::<PcurveHolder>::get_owned(&table, 12).unwrap();
    let curve = PCurve::try_from(&source).unwrap();
    let (u, v) = (0.4_f64, 0.2_f64);
    let radius = 1.55 + 5.0 * v.cos();
    assert!(
        curve.subs(0.0).distance(Point3::new(
            radius * u.cos(),
            radius * u.sin(),
            5.0 * v.sin()
        )) < 1e-12
    );
}

#[test]
fn signed_torus_reverses_bounds_and_preserves_declared_face_sense() {
    use truck_stepio::r#in::step_geometry::*;
    let table =
        Table::from_step(include_str!("input/fixtures/negative-torus-synthetic.step")).unwrap();
    let (shell, losses) = table
        .to_compressed_shell_with_losses(102, &table.shell[&102])
        .unwrap();
    assert!(losses.is_empty(), "{losses:?}");
    assert_eq!(shell.faces.len(), 2);
    let positive = Table::from_step(
        &include_str!("input/fixtures/negative-torus-synthetic.step").replace("-1.55", "1.55"),
    )
    .unwrap();
    let original = positive
        .to_compressed_shell(102, &positive.shell[&102])
        .unwrap();
    for (signed, unsigned) in shell.faces.iter().zip(&original.faces) {
        let expected: Vec<_> = unsigned.boundaries[0]
            .iter()
            .rev()
            .map(|edge| (edge.index, !edge.orientation))
            .collect();
        let actual: Vec<_> = signed.boundaries[0]
            .iter()
            .map(|edge| (edge.index, edge.orientation))
            .collect();
        assert_eq!(actual, expected);
    }
    let u = 0.4_f64;
    let v = 0.2_f64;
    let radius = 1.55 + 5.0 * v.cos();
    let expected = Point3::new(radius * u.cos(), radius * u.sin(), 5.0 * v.sin());
    let normal = Vector3::new(v.cos() * u.cos(), v.cos() * u.sin(), v.sin());
    assert!(shell.faces[0].surface.subs(u, v).distance(expected) < 1e-12);
    assert!((shell.faces[0].surface.normal(u, v) - normal).magnitude() < 1e-12);
    assert!(shell.faces[1].surface.subs(v, u).distance(expected) < 1e-12);
    assert!((shell.faces[1].surface.normal(v, u) + normal).magnitude() < 1e-12);
}
