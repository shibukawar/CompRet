# route_extraction

Extract, format, and rank reaction routes given a proof tree.

## Requirments

- networkx
- rdkit
- tqdm

## Behavior

- Extract routes with "result_dir/pt.txt", "result_dir/finalResult.txt". result_dir is scpecified with "setting.txt"
- Create new directory "route" under the result_dir.
- Convert each route into two files, "state.txt" and "reaction.txt".
- Create folders each contains the above two files and add them under "(result_dir)/route/"
- Output ranking result to "log.txt"
