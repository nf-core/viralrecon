import xml.etree.ElementTree as ET
import pandas as pd
import argparse
import json
import re
import sys


# ---------------------- Argument parsing ----------------------

def parser_args(args=None):
    parser = argparse.ArgumentParser(
        description="Create JSON file from XML with HIV mutation information.",
        epilog="Example: python xml_processing.py -x HIVDB_9.8.xml -o comments.json"
    )
    parser.add_argument(
        "-x", "--xml_file", type=str, default="./HIVDB_9.8.xml",
        help="Path to input XML file (default: ./HIVDB_9.8.xml)."
    )
    parser.add_argument(
        "-oj", "--output_json", type=str, default="./comments.json",
        help="Path to output JSON file (default: ./comments.json)."
    )
    parser.add_argument(
        "-oc", "--output_csv", type=str, default="./comments.csv",
        help="Path to output CSV file."
    )
    return parser.parse_args(args)


# ---------------------- Load XML (including comments) ----------------------

def load_xml_root_with_comments(xml_path):
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    tree = ET.parse(xml_path, parser)
    return tree.getroot()


# ---------------------- Gene / Drug-class mapping ----------------------

def extract_gene_drugclass(root):
    gene_drugclass = {}

    for gd in root.findall(".//GENE_DEFINITION"):
        gene = gd.findtext("NAME")
        drug_classes = gd.findtext("DRUGCLASSLIST")
        if drug_classes:
            gene_drugclass[gene] = [dc.strip() for dc in drug_classes.split(",")]

    return gene_drugclass


def make_reverse_drugclass_map(gene_drugclass):
    """Create mapping drugclass → gene."""
    reverse = {}
    for gene, classes in gene_drugclass.items():
        for dc in classes:
            reverse[dc] = gene
    return reverse


# ---------------------- Mutation parsing ----------------------

def expand_mutation(mut_str):
    """Expand mutation strings like E138K/A/T → [E138K, E138A, E138T]."""

    m = re.match(r"^([A-Z])([0-9]+)([A-Z/]+)$", mut_str)
    if not m:
        return [mut_str]

    ref_aa, pos, alts = m.groups()
    alt_list = alts.split("/") if "/" in alts else list(alts)

    return [f"{ref_aa}{pos}{aa}" for aa in alt_list]


# ---------------------- Extract comments and mutations ----------------------

def extract_mutations_from_comments(root, gene_drugclass):
    MUTATION_REGEX = re.compile(r"^([A-Z][0-9]+[A-Z/]+)")
    DRUGCLASS_REGEX = re.compile(r"<DRUGCLASS>(.*?)</DRUGCLASS>")

    comments = {}
    drugclass_to_gene = make_reverse_drugclass_map(gene_drugclass)

    for c in root.findall(".//COMMENT_STRING"):
        cid = c.get("id")

        text_el = c.find("TEXT")
        if text_el is None:
            continue

        comment_text = text_el.text.strip()

        # Extract DRUGCLASS (inside XML comments)
        raw_xml = ET.tostring(c, encoding="unicode")
        m_dc = DRUGCLASS_REGEX.search(raw_xml)
        drug_class = m_dc.group(1).strip() if m_dc else None
        gene = drugclass_to_gene.get(drug_class)

        # Mutation detection
        m_mut = MUTATION_REGEX.match(comment_text)
        if m_mut:
            mut_list = expand_mutation(m_mut.group(1))
        else:
            mut_list = [f"{drug_class} comments" if drug_class else "Unassigned comments"]

        comments[cid] = {
            "gene": gene,
            "drug_class": drug_class,
            "comment": comment_text,
            "mutation": mut_list
        }

    return comments


# ---------------------- Reorganize by mutation key ----------------------

def reorganize_by_mutation(comments_by_id):
    mutations = {}

    for cid, entry in comments_by_id.items():
        gene = entry["gene"]
        drug_class = entry["drug_class"]
        comment = entry["comment"]

        for mut in entry["mutation"]:
            if mut not in mutations:
                mutations[mut] = {
                    "gene": gene,
                    "drug_class": drug_class,
                    "comments": {cid: comment}
                }
            else:
                # Avoid duplicates
                existing = mutations[mut]["comments"]
                if cid not in existing and comment not in existing.values():
                    existing[cid] = comment

    return mutations


# ---------------------- Output ----------------------

def save_json(data, output_path):
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


# ---------------------- JSON to table ----------------------

def json_to_table(json_file):
    rows = []
    for mutation, info in json_file.items():
        gene = info.get("gene")
        drug_class = info.get("drug_class")
        comments = info.get("comments", {})

        for comment_id, comment_text in comments.items():
            rows.append({
                "Mutation/Genes": mutation,
                "Gene": gene,
                "Drug Class": drug_class,
                "Comment ID": comment_id,
                "Comment": comment_text
            })

    df = pd.DataFrame(rows)
    return df

# ---------------------- Main ----------------------

def main(args=None):
    args = parser_args(args)

    root = load_xml_root_with_comments(args.xml_file)
    gene_drugclass = extract_gene_drugclass(root)
    comments = extract_mutations_from_comments(root, gene_drugclass)
    mutations = reorganize_by_mutation(comments)

    save_json(mutations, args.output_json)
    print(f"JSON file generated: {args.output_json}")

    df = json_to_table(mutations)
    df.to_csv(args.output_csv, index=False, encoding="utf-8")
    print(f"CSV file generated: {args.output_csv}")

if __name__ == "__main__":
    sys.exit(main())
