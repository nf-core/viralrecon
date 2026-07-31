# Plan: Fix Nextflow Strict Syntax Lint Errors in viralrecon

## Context

Nextflow strict syntax compliance is becoming mandatory for nf-core pipelines (spring 2026 deadline). `nextflow lint` on viralrecon produces **49 errors** and **99 warnings**. We're fixing **local files only** (not `modules/nf-core/*` or `subworkflows/nf-core/*`).

That leaves **31 errors** and ~75 warnings to fix across local files.

---

## Phase 1: Trivial Fixes (Completed)

### 1a. Missing config file (1 error)

- **File:** `nextflow.config:263`
- **Error:** `Invalid include source: conf/test_full_sispa.config` — file doesn't exist
- **Fix:** Create an empty/stub `conf/test_full_sispa.config` (or copy from `test_sispa.config` as a base), or remove the `test_full_sispa` profile entry

### 1b. Variable shadowing in `modules/local/sierralocal/main.nf` (6 errors)

- **Lines 28-33:** `def hivdb_xml`, `def apobec_drm`, etc. shadow `input:` parameter names
- **Fix:** Rename script-block vars → `hivdb_xml_arg`, `apobec_drm_arg`, `apobec_csv_arg`, `unusual_csv_arg`, `sdrms_csv_arg`, `mutation_csv_arg`; update references in the script string below

### 1c. Variable shadowing in `subworkflows/local/consensus_bcftools/main.nf` (3 errors)

- **Line 73:** `.map { meta, vcf, tbi, fasta -> ... }` shadows `take:` parameter names
- **Fix:** Rename closure params → `.map { meta, vcf_file, tbi_file, fasta_file -> ... }`

### 1d. `Nextflow.error()` → `error()` in `main.nf` (3 errors)

- **Lines 176, 186, 195:** `Nextflow.error(...)` not valid in strict mode
- **Fix:** Replace with `error(...)` (built-in function)

---

## Phase 2: Moderate Fixes

### 2a. `main.nf` top-level statements (7 errors + downstream errors)

- **Lines 18-41:** `def`, `if`, `params.*=` mixed with `include`/`workflow` declarations
- **Downstream:** 6 more errors (`VIRALRECON`/`artic_scheme`/`primer_set` not defined) resolve once the parsing succeeds
- **Fix:**
  1. Move lines 18-19 (`def primer_set`, `def primer_set_version`) and lines 22-41 (validation + params assignment) into the `getGenomeAttribute` area or a new `initGenomeParams()` function
  2. Call that function from within the entry `workflow {}` block, storing results
  3. Pass `artic_scheme` through to `NFCORE_VIRALRECON` via the workflow call
  4. Move `params.artic_scheme = ...` (line 38) to config or handle without mutating params (W5 warning too)
  5. Move `params.fasta = ...` etc. (lines 43-50) into config or a function

### 2b. Conditional `include` in `workflows/viralrecon.nf` (2 errors)

- **Line 91-95:** `if (params.platform == 'illumina') { include { PREPARE_GENOME_ILLUMINA ... } }`
- **Line 92:** `Unexpected input: 'include'`
- **Fix:** Move both includes to top-level (unconditional):
  ```groovy
  include { PREPARE_GENOME_ILLUMINA } from '../subworkflows/local/prepare_genome_illumina'
  include { PREPARE_GENOME_NANOPORE } from '../subworkflows/local/prepare_genome_nanopore'
  ```
  Then branch at the call site:
  ```groovy
  if (params.platform == 'illumina') {
      PREPARE_GENOME_ILLUMINA(...)
  } else {
      PREPARE_GENOME_NANOPORE(...)
  }
  ```
  Note: can no longer alias both as `PREPARE_GENOME` — downstream references need updating.

### 2c. `WorkflowCommons` not defined (2 errors)

- **Files:** `subworkflows/local/variants_ivar/main.nf:43`, `subworkflows/local/variants_bcftools/main.nf:40`
- **Problem:** `lib/` Groovy classes aren't auto-imported in strict mode
- **Fix:** Convert the two referenced static methods (`getNumLinesInFile`, `getNumVariantsFromBCFToolsStats`) to regular Nextflow functions defined in the subworkflow files (or a shared `.nf` file that's `include`'d). The functions are simple enough to inline:
  ```groovy
  def getNumLinesInFile(input_file) {
      def num_lines = 0
      input_file.eachLine { line -> num_lines++ }
      return num_lines
  }
  ```
  Note: `workflows/viralrecon.nf` has many more `WorkflowCommons` references (~19), but the lint only flags the 2 in subworkflows. The viralrecon.nf ones may be masked by the parse failure (error #4). Once #4 is fixed, more `WorkflowCommons` errors will likely appear — all will need the same treatment (convert to Nextflow functions).

### 2d. `new JsonSlurper()` in `subworkflows/local/fastq_trim_fastp_fastqc/main.nf` (1 error)

- **Line 15:** `new JsonSlurper()` — `new` keyword disallowed in strict syntax
- **Fix:** Replace with Nextflow's built-in JSON parsing. Use `jsonSlurper()` factory or `file.text` with `groovy.json.JsonSlurper` imported differently. Alternatively, move the JSON parsing into the process script block (bash `jq`).

---

## Phase 3: Config File Restructuring (18 errors)

- **Files:** `conf/modules_illumina.config`, `conf/modules_nanopore.config`, `conf/modules.config`
- **Problem:** `if`/variable statements mixed with config declarations
- **Fix:** Requires significant restructuring of config files

---

## Phase 4: Warnings (local files only, ~75 warnings)

### 4a. Implicit `it` → explicit closure parameter (~35 warnings in local files)

- Across ~14 local subworkflow files
- Mechanical: `{ it[1] }` → `{ v -> v[1] }`, `{ it }` → `{ v -> v }`, etc.
- Files: `additional_annotation`, `assembly_qc`, `consensus_qc`, `filter_bam_samtools`, `hiv_resitance_detection`, `prepare_genome_illumina`, `prepare_genome_nanopore`, `variants_bcftools`, `variants_long_table`, `utils_nfcore_viralrecon_pipeline`

### 4b. Unused parameters → prefix with `_` (~23 warnings in local files)

- Across ~8 local subworkflow files
- Mechanical: `meta` → `_meta`, `illumina` → `_illumina`, etc. in `.filter`/`.map` closures
- Files: `assembly_minia`, `assembly_spades`, `assembly_unicycler`, `variants_bcftools`, `variants_ivar`, `variants_qc`, `utils_nfcore_viralrecon_pipeline`

### 4c. Unused variables (~2 warnings in local files)

- `main.nf:41` — `def artic_scheme` (will be addressed in Phase 2a)
- `utils_nfcore_viralrecon_pipeline/main.nf:297` — declared but unused variable

---

## Implementation Order

1. ~~Phase 1a-1d (trivial fixes) — 13 errors resolved~~
2. Phase 2b (conditional includes in viralrecon.nf) — 2 errors resolved
3. Phase 2a (main.nf restructuring) — 13 errors resolved (7 direct + 6 downstream)
4. Phase 2c (WorkflowCommons) — 2 errors resolved (+ likely more uncovered)
5. Phase 2d (JsonSlurper) — 1 error resolved
6. Phase 3 (config restructuring) — 18 errors resolved
7. Phase 4a-4c (warnings) — ~60 warnings resolved

## Verification

After each phase, run:

```bash
nextflow lint . -exclude ".git,.nf-test,nf-test.config,work"
```
