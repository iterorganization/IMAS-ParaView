def test_validate_dummy_ids(dummy_ids):
    """Validates the dummy IDS object created by the fixture."""

    ids = dummy_ids
    ids.validate()
