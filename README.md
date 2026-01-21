# neutropenia_redcap
Ye olde specious harcoded example
```
uv run python -m neutropenia_redcap.convert \
    --data_location ~/r/Eli/redcap_resources/post_processed_full_dataset.tsv \
    --output_dir . \
    --input_format RAW_TSV \
    --output_format REDCAP
```
