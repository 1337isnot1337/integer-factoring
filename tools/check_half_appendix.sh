#!/usr/bin/env bash
set -euo pipefail

# Build the pinned Appendix B development independently of the factoring
# reference: the two upstream projects currently use different Lean versions.
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="$root/.lake/appendix-source"
revision="16a5ce41304fcdaea4e0be7c132e1c361dc26da5"
logs="$root/build/half-appendix"
mkdir -p "$root/.lake" "$logs"
if [[ ! -d "$source_dir" ]]; then
  git clone --filter=blob:none --no-checkout https://github.com/anthropics/formal-math.git "$source_dir"
  git -C "$source_dir" checkout --detach "$revision"
fi
if [[ "$(git -C "$source_dir" rev-parse HEAD)" != "$revision" ]] ||
   [[ -n "$(git -C "$source_dir" status --porcelain --untracked-files=no)" ]]; then
  echo "Appendix source is not the clean pinned revision" >&2
  exit 1
fi
cd "$source_dir/eth"
lean --version | tee "$logs/toolchain.log"
lake exe cache get > "$logs/cache.log" 2>&1
lake build GeneralPieces.Main CountKSat.Host > "$logs/source-build.log" 2>&1

# The stage is below the source project so Lean can set a separate root while
# still finding the source project's built imports through `lake env`.
stage="$(mktemp -d "$source_dir/eth/.lake/half-bridge.XXXXXXXX")"
reference="$root/.lake/upstream-reference"
reference_revision="c91cc3f4947caca8977e81a5d14244f73b5ae2b1"
if [[ ! -d "$reference" ]]; then
  git clone --filter=blob:none --no-checkout https://github.com/IsaacRudich/SubexponentialDeterministicFactoring.git "$reference"
  git -C "$reference" checkout --detach "$reference_revision"
fi
if [[ "$(git -C "$reference" rev-parse HEAD)" != "$reference_revision" ]] ||
   [[ -n "$(git -C "$reference" status --porcelain --untracked-files=no)" ]]; then
  echo "Factoring source is not the clean pinned revision" >&2
  exit 1
fi
port="$stage/factoring"
mkdir -p "$port/.lake"
git -C "$reference" archive "$reference_revision" | tar -x -C "$port"
cp "$source_dir/eth/lean-toolchain" "$port/lean-toolchain"
cp "$source_dir/eth/lake-manifest.json" "$port/lake-manifest.json"
ln -s "$source_dir/eth/.lake/packages" "$port/.lake/packages"
cat > "$port/lakefile.toml" <<'EOF'
name = "factoring_port"

[[require]]
name = "mathlib"
git = "https://github.com/leanprover-community/mathlib4"
rev = "0df444a360eaa60ab8c11dca51a86af692955474"

[[lean_lib]]
name = "FactoringAux"
globs = ["FactoringAux.+"]

[[lean_lib]]
name = "Challenge"
EOF
(cd "$port" && lake build FactoringAux.MulCount FactoringAux.Cluster FactoringAux.EngineGroup FactoringAux.EngineCostRate FactoringAux.Search Challenge) > "$logs/factoring-build.log" 2>&1
cp "$root/formal/HalfPad.lean" "$stage/HalfPad.lean"
cp "$root/formal-appendix/ExplicitPad.lean" "$stage/ExplicitPad.lean"
cp "$root/formal-appendix/Tabulate.lean" "$stage/Tabulate.lean"
cp "$root/formal-appendix/CountBridge.lean" "$stage/CountBridge.lean"
cp "$root/formal-appendix/CliqueNetwork.lean" "$stage/CliqueNetwork.lean"
cp "$root/formal-appendix/CliqueRows.lean" "$stage/CliqueRows.lean"
cp "$root/formal-appendix/RowSimplify.lean" "$stage/RowSimplify.lean"
cp "$root/formal-appendix/TopValue.lean" "$stage/TopValue.lean"
cp "$root/formal-appendix/TopEvaluator.lean" "$stage/TopEvaluator.lean"
cp "$root/formal-appendix/TopGrouping.lean" "$stage/TopGrouping.lean"
cp "$root/formal-appendix/SelectedTop.lean" "$stage/SelectedTop.lean"
cp "$root/formal-appendix/FastTop.lean" "$stage/FastTop.lean"
cp "$root/formal-appendix/TopArithmetic.lean" "$stage/TopArithmetic.lean"
cp "$root/formal-appendix/TopStates.lean" "$stage/TopStates.lean"
cp "$root/formal-appendix/TopStateArithmetic.lean" "$stage/TopStateArithmetic.lean"
cp "$root/formal-appendix/EntryBound.lean" "$stage/EntryBound.lean"
cp "$root/formal-appendix/BinaryFac.lean" "$stage/BinaryFac.lean"
cp "$root/formal-appendix/BinaryNet.lean" "$stage/BinaryNet.lean"
cp "$root/formal-appendix/BinaryCoarsen.lean" "$stage/BinaryCoarsen.lean"
cp "$root/formal-appendix/SelectedEntry.lean" "$stage/SelectedEntry.lean"
cp "$root/formal-appendix/BoundedTop.lean" "$stage/BoundedTop.lean"
cp "$root/formal-appendix/HalfSchedule.lean" "$stage/HalfSchedule.lean"
cp "$root/formal-appendix/PadCost.lean" "$stage/PadCost.lean"
cp "$root/formal-appendix/SelectedPadCost.lean" "$stage/SelectedPadCost.lean"
cp "$root/formal/HalfPrefix.lean" "$stage/HalfPrefix.lean"
cp "$root/formal-appendix/PrefixEvaluator.lean" "$stage/PrefixEvaluator.lean"
cp "$root/formal-appendix/PureFactor.lean" "$stage/PureFactor.lean"
cp "$root/formal-appendix/LevelSelect.lean" "$stage/LevelSelect.lean"
cp "$root/formal-appendix/SelectedRate.lean" "$stage/SelectedRate.lean"
cp "$root/formal-appendix/WrittenPoly.lean" "$stage/WrittenPoly.lean"
cp "$root/formal-appendix/InputSize.lean" "$stage/InputSize.lean"
cp "$root/formal-appendix/FragmentCost.lean" "$stage/FragmentCost.lean"
cp "$root/formal-appendix/EvaluatorSteps.lean" "$stage/EvaluatorSteps.lean"
cp "$root/formal-appendix/BoundedLength.lean" "$stage/BoundedLength.lean"
cp "$root/formal-appendix/PrefixSteps.lean" "$stage/PrefixSteps.lean"
cp "$root/formal-appendix/EvaluatorRate.lean" "$stage/EvaluatorRate.lean"
cp "$root/formal-appendix/ReferenceStatement.lean" "$stage/ReferenceStatement.lean"
cp "$root/formal-appendix/CertifiedPrefix.lean" "$stage/CertifiedPrefix.lean"
cp "$root/formal-appendix/SearchPoly.lean" "$stage/SearchPoly.lean"
cp "$root/formal-appendix/PaddedWordMachine.lean" "$stage/PaddedWordMachine.lean"
cp "$root/formal-appendix/GroupingSearchPoly.lean" "$stage/GroupingSearchPoly.lean"
cp "$root/formal-appendix/GroupingSearchCost.lean" "$stage/GroupingSearchCost.lean"
cp "$root/formal-appendix/LevelPoly.lean" "$stage/LevelPoly.lean"
cp "$root/formal-appendix/PaddedMachineRate.lean" "$stage/PaddedMachineRate.lean"
cp "$root/formal-appendix/SizeParamPoly.lean" "$stage/SizeParamPoly.lean"
cp "$root/formal-appendix/TopParamsPoly.lean" "$stage/TopParamsPoly.lean"
cp "$root/formal-appendix/MasterPad.lean" "$stage/MasterPad.lean"
cp "$root/formal-appendix/GroupingCodePoly.lean" "$stage/GroupingCodePoly.lean"
cp "$root/formal-appendix/CrossQueryPoly.lean" "$stage/CrossQueryPoly.lean"
cp "$root/formal-appendix/CrossQueryCorrect.lean" "$stage/CrossQueryCorrect.lean"
cp "$root/formal-appendix/GroupingPredicate.lean" "$stage/GroupingPredicate.lean"
lake env bash -c '
  set -euo pipefail
  export LEAN_PATH="$1:$1/factoring/.lake/build/lib/lean:$LEAN_PATH"
  lean --root="$1" -o "$1/HalfPad.olean" "$1/HalfPad.lean"
  lean --root="$1" -o "$1/ExplicitPad.olean" "$1/ExplicitPad.lean"
  lean --root="$1" -o "$1/Tabulate.olean" "$1/Tabulate.lean"
  lean --root="$1" "$1/CountBridge.lean"
  lean --root="$1" -o "$1/CliqueNetwork.olean" "$1/CliqueNetwork.lean"
  lean --root="$1" -o "$1/CliqueRows.olean" "$1/CliqueRows.lean"
  lean --root="$1" -o "$1/RowSimplify.olean" "$1/RowSimplify.lean"
  lean --root="$1" -o "$1/TopValue.olean" "$1/TopValue.lean"
  lean --root="$1" -o "$1/TopEvaluator.olean" "$1/TopEvaluator.lean"
  lean --root="$1" -o "$1/TopGrouping.olean" "$1/TopGrouping.lean"
  lean --root="$1" -o "$1/SelectedTop.olean" "$1/SelectedTop.lean"
  lean --root="$1" -o "$1/FastTop.olean" "$1/FastTop.lean"
  lean --root="$1" -o "$1/TopArithmetic.olean" "$1/TopArithmetic.lean"
  lean --root="$1" -o "$1/TopStates.olean" "$1/TopStates.lean"
  lean --root="$1" -o "$1/TopStateArithmetic.olean" "$1/TopStateArithmetic.lean"
  lean --root="$1" -o "$1/EntryBound.olean" "$1/EntryBound.lean"
  lean --root="$1" -o "$1/BinaryFac.olean" "$1/BinaryFac.lean"
  lean --root="$1" -o "$1/BinaryNet.olean" "$1/BinaryNet.lean"
  lean --root="$1" -o "$1/BinaryCoarsen.olean" "$1/BinaryCoarsen.lean"
  lean --root="$1" -o "$1/SelectedEntry.olean" "$1/SelectedEntry.lean"
  lean --root="$1" -o "$1/BoundedTop.olean" "$1/BoundedTop.lean"
  lean --root="$1" -o "$1/HalfSchedule.olean" "$1/HalfSchedule.lean"
  lean --root="$1" -o "$1/PadCost.olean" "$1/PadCost.lean"
  lean --root="$1" -o "$1/SelectedPadCost.olean" "$1/SelectedPadCost.lean"
  lean --root="$1" -o "$1/HalfPrefix.olean" "$1/HalfPrefix.lean"
  lean --root="$1" -o "$1/PrefixEvaluator.olean" "$1/PrefixEvaluator.lean"
  lean --root="$1" -o "$1/PureFactor.olean" "$1/PureFactor.lean"
  lean --root="$1" -o "$1/LevelSelect.olean" "$1/LevelSelect.lean"
  lean --root="$1" "$1/SelectedRate.lean"
  lean --root="$1" -o "$1/WrittenPoly.olean" "$1/WrittenPoly.lean"
  lean --root="$1" -o "$1/InputSize.olean" "$1/InputSize.lean"
  lean --root="$1" -o "$1/FragmentCost.olean" "$1/FragmentCost.lean"
  lean --root="$1" -o "$1/EvaluatorSteps.olean" "$1/EvaluatorSteps.lean"
  lean --root="$1" -o "$1/BoundedLength.olean" "$1/BoundedLength.lean"
  lean --root="$1" -o "$1/PrefixSteps.olean" "$1/PrefixSteps.lean"
  lean --root="$1" -o "$1/EvaluatorRate.olean" "$1/EvaluatorRate.lean"
  lean --root="$1" "$1/ReferenceStatement.lean"
  lean --root="$1" "$1/CertifiedPrefix.lean"
  lean --root="$1" "$1/SearchPoly.lean"
  lean --root="$1" -o "$1/PaddedWordMachine.olean" "$1/PaddedWordMachine.lean"
  lean --root="$1" "$1/GroupingSearchPoly.lean"
  lean --root="$1" -o "$1/GroupingSearchCost.olean" "$1/GroupingSearchCost.lean"
  lean --root="$1" -o "$1/LevelPoly.olean" "$1/LevelPoly.lean"
  lean --root="$1" "$1/PaddedMachineRate.lean"
  lean --root="$1" -o "$1/SizeParamPoly.olean" "$1/SizeParamPoly.lean"
  lean --root="$1" -o "$1/TopParamsPoly.olean" "$1/TopParamsPoly.lean"
  lean --root="$1" "$1/MasterPad.lean"
  lean --root="$1" -o "$1/GroupingCodePoly.olean" "$1/GroupingCodePoly.lean"
  lean --root="$1" -o "$1/CrossQueryPoly.olean" "$1/CrossQueryPoly.lean"
  lean --root="$1" -o "$1/CrossQueryCorrect.olean" "$1/CrossQueryCorrect.lean"
  lean --root="$1" "$1/GroupingPredicate.lean"
' _ "$stage" > "$logs/bridge.log" 2>&1
python3 - "$logs/bridge.log" <<'PY'
import re
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text()
reports = re.findall(r"'([^']+)' depends on axioms:\s*\[([^\]]*)\]", text, re.S)
expected = "IntegerFactoring.HalfAppendix.exists_explicit_padded_evaluator"
count_bridge = "IntegerFactoring.HalfAppendix.evaluator_counts_of_value_bridge"
table = "IntegerFactoring.HalfAppendix.evaluator_counts_tabulated_of_value_bridge"
clique = "IntegerFactoring.HalfAppendix.cliqueNetwork_value"
top = "IntegerFactoring.HalfAppendix.evaluator_counts_top"
selected = "IntegerFactoring.HalfAppendix.evaluator_counts_selected_top"
row = "IntegerFactoring.HalfAppendix.cliqueQuery_bound"
fast = "IntegerFactoring.HalfAppendix.evaluator_counts_selected_fast"
width = "IntegerFactoring.HalfAppendix.selected_width_le"
states = "IntegerFactoring.HalfAppendix.selectedTop_state_card"
state_rate = "IntegerFactoring.HalfAppendix.selected_state_exponent_le"
entry = "IntegerFactoring.HalfAppendix.clique_entry_le_card"
explicit = "IntegerFactoring.HalfAppendix.evaluator_counts_explicit_top"
pad_rate = "IntegerFactoring.HalfAppendix.scheduled_pad_sqrt_cost"
pure_factor = "IntegerFactoring.HalfAppendix.pure_factor_all"
level_rate = "IntegerFactoring.HalfAppendix.selected_level_pad_sqrt_cost"
written_poly = "IntegerFactoring.HalfAppendix.written_evaluator_polyT"
prefix_steps = "IntegerFactoring.HalfAppendix.prefix_evaluator_correct_steps"
reference_spec = "IntegerFactoring.HalfAppendix.pure_factor_reference_spec"
certified_prefix = "IntegerFactoring.HalfAppendix.prefix_evaluator_sqrt_steps"
search_poly = "IntegerFactoring.HalfAppendix.search_polyT"
machine = "IntegerFactoring.HalfAppendix.padded_factoring_machine"
group_search = "IntegerFactoring.HalfAppendix.padded_firstGood_of_count_polyT"
group_cost = "IntegerFactoring.HalfAppendix.scheduled_grouping_pad_cost"
level_poly = "IntegerFactoring.HalfAppendix.selectedLevel_polyT"
machine_rate = "IntegerFactoring.HalfAppendix.padded_factoring_exp"
size_poly = "IntegerFactoring.HalfAppendix.sOf_unary"
params_poly = "IntegerFactoring.HalfAppendix.selected_padExponent_unary"
master_cost = "IntegerFactoring.HalfAppendix.masterExponent_sqrt_cost"
group_sizes = "IntegerFactoring.HalfAppendix.groupSizesOK_eq"
group_predicate = "IntegerFactoring.HalfAppendix.goodCode_polyT"
group_correct = "IntegerFactoring.HalfAppendix.goodCode_eq_okG"
if "error:" in text or "sorryAx" in text or not {expected, count_bridge, table, clique, top, selected, row, fast, width, states, state_rate, entry, explicit, pad_rate, pure_factor, level_rate, written_poly, prefix_steps, reference_spec, certified_prefix, search_poly, machine, group_search, group_cost, level_poly, machine_rate, size_poly, params_poly, master_cost, group_sizes, group_predicate, group_correct} <= {name for name, _ in reports}:
    raise SystemExit("Missing or failed Appendix bridge axiom reports")
allowed = {"propext", "Classical.choice", "Quot.sound"}
for name, axioms in reports:
    actual = {a.strip() for a in axioms.split(",") if a.strip()}
    if actual - allowed:
        raise SystemExit(f"Unexpected axioms for {name}: {sorted(actual)}")
print(f"Checked {len(reports)} axiom reports; logs: {sys.argv[1]}")
PY
