from chemparse import parse_formula

input_file = "composition_list"
output_file = "chemparse_result"

parsed_compositions = []

with open(input_file, "r") as f:
    for line in f:
        formula = line.strip()
        if not formula:
            continue
        try:
            parsed = parse_formula(formula)
            parsed_compositions.append((formula, parsed))
        except Exception as e:
            parsed_compositions.append((formula, f"Error: {e}"))

with open(output_file, "w") as f:
    for formula, result in parsed_compositions:
        f.write(f"{formula}: {result}\n")

print(f"{len(parsed_compositions)} formulas parsed and saved to '{output_file}'")

