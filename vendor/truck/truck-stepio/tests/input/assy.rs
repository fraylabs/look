use truck_stepio::r#in::*;

const STEP_DIRECTORY: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../resources/step/");

#[test]
fn occt_assy() {
    let step_string = std::fs::read_to_string([STEP_DIRECTORY, "occt-assy.step"].concat()).unwrap();
    let table = Table::from_step(&step_string).unwrap();
    let assy = table.step_assy().unwrap();
    assert_eq!(assy.len(), 3);
    let top = assy.top_nodes().next().unwrap();
    let paths = assy.maximal_paths_iter(top.index()).collect::<Vec<_>>();
    assert_eq!(paths.len(), 5);
}

#[test]
fn nullable_assembly_metadata_preserves_both_occurrences() {
    let table = Table::from_step(include_str!("fixtures/assembly-pds-dollar.step")).unwrap();
    assert_eq!(table.next_assembly_usage_occurrence.len(), 2);
    assert_eq!(table.product_definition_shape.len(), 4);
    assert_eq!(table.item_defined_transformation.len(), 2);
    let assy = table
        .step_assy()
        .expect("nullable metadata must retain the graph");
    assert_eq!(assy.all_edges().count(), 2);
    assert_eq!(assy.top_nodes().count(), 1);
}
