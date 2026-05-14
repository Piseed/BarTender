import argparse
import json
import os
import sys
import edlib
import pandas as pd
from Bio import SeqIO

# --- Прием параметров от CLI менеджера ---
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_fastq", required=True)
parser.add_argument("-m", "--manifest", required=True)
parser.add_argument("-d", "--dbdir", required=True)
parser.add_argument("-o", "--outdir", required=True)
parser.add_argument("--max_errors_adapter", type=int, required=True)
parser.add_argument("--max_errors_barcode", type=int, required=True)
parser.add_argument("--indel_buffer", type=int, required=True)
args = parser.parse_args()

MAX_ERRORS_ADAPTER = args.max_errors_adapter
MAX_ERRORS_BARCODE = args.max_errors_barcode
INDEL_BUFFER = args.indel_buffer

# --- Парсинг адаптеров ---
json_path = os.path.join(args.dbdir, "info_adapters.json")
if not os.path.exists(json_path):
    sys.exit(f"Ошибка: Файл '{json_path}' не найден.")

with open(json_path, "r") as jf:
    adapter_config = json.load(jf)

a5_seq = adapter_config.get("adapterf", "").upper()
a3_seq = adapter_config.get("adapterr", "").upper()
IS_DUAL_INDEX = bool(a5_seq and a3_seq)

# Настройка пулов поиска в зависимости от типа библиотеки
ADAPTERS_5PRIME = {}
if IS_DUAL_INDEX:
    ADAPTERS_5PRIME["F_raw"] = a5_seq
    ADAPTERS_5PRIME["R_raw"] = a3_seq
else:
    # Для Single-index на 5'-конце рида мы ищем ТОЛЬКО прямой адаптер (Forward нить)
    ADAPTERS_5PRIME["F_raw"] = a5_seq

PAIRED_3PRIME = {}
if IS_DUAL_INDEX:
    PAIRED_3PRIME = {
        "F_raw": ("R_revcomp", a3_seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]),
        "R_raw": ("F_revcomp", a5_seq.translate(str.maketrans("ATCG", "TAGC"))[::-1])
    }

# --- 1. Загрузка коллекций баркодов ---
barcode_db = {}
groups = ["F_raw", "F_revcomp"] if not IS_DUAL_INDEX else ["F_raw", "F_revcomp", "R_raw", "R_revcomp"]

for group in groups:
    barcode_db[group] = {}
    fasta_name = os.path.join(args.dbdir, f"collection_{group}.fasta")
    if os.path.exists(fasta_name):
        for r in SeqIO.parse(fasta_name, "fasta"):
            clean_id = r.id.rsplit('_', 1)[0] if '_' in r.id else r.id
            barcode_db[group][str(r.seq).upper()] = clean_id

# --- 2. Загрузка Excel манифеста ---
df = pd.read_excel(args.manifest)
sample_manifest = {}
for _, row in df.iterrows():
    b1 = str(row['Barcode1']).strip()
    if IS_DUAL_INDEX and 'Barcode2' in df.columns and pd.notna(row['Barcode2']):
        key = frozenset((b1, str(row['Barcode2']).strip()))
    else:
        key = b1
    sample_manifest[key] = str(row['SampleID']).strip()

# --- 3. Инициализация файлов вывода ---
os.makedirs(os.path.join(args.outdir, "demux_fastq"), exist_ok=True)
out_files = {"unidentified": open(os.path.join(args.outdir, "demux_fastq/unidentified.fastq"), "w")}
report = open(os.path.join(args.outdir, "demux_report.tsv"), "w")

if IS_DUAL_INDEX:
    report.write("read_id\tstatus\tadapter_5p\tbc1_id\tadapter_3p\tbc2_id\tsample_id\n")
else:
    report.write("read_id\tstatus\tfound_adapter\torientation\tbc_id\tsample_id\n")

count, fully_identified = 0, 0

# --- 4. Основной пайплайн обработки ---
for record in SeqIO.parse(args.input_fastq, "fastq"):
    read_seq = str(record.seq).upper()
    read_id = record.id
    count += 1
    
    if IS_DUAL_INDEX:
        # =====================================================================
        # ЛОГИКА ДЛЯ DUAL-INDEX (Остается неизменной)
        # =====================================================================
        best_5p_name, best_5p_err = None, MAX_ERRORS_ADAPTER + 1
        best_5p_start, best_5p_end = 0, 0
        
        half_read = read_seq[:len(read_seq)//2]
        for name, seq in ADAPTERS_5PRIME.items():
            res = edlib.align(seq, half_read, mode="HW", task="locations")
            if res["editDistance"] != -1 and res["editDistance"] <= MAX_ERRORS_ADAPTER:
                if res["editDistance"] < best_5p_err:
                    best_5p_err, best_5p_name = res["editDistance"], name
                    best_5p_start, best_5p_end = res["locations"][0]
                    
        if not best_5p_name:
            SeqIO.write(record, out_files["unidentified"], "fastq")
            report.write(f"{read_id}\tUnidentifiable_5p\tNA\tNA\tNA\tNA\tNA\n")
            continue

        bc1_id, bc1_err = None, MAX_ERRORS_BARCODE + 1
        window_bc1_end = min(len(read_seq), best_5p_start + INDEL_BUFFER)
        sub_read_bc1 = read_seq[0 : window_bc1_end]
        for bc_seq, bc_id_candidate in barcode_db[best_5p_name].items():
            res_bc = edlib.align(bc_seq, sub_read_bc1, mode="HW")
            if res_bc["editDistance"] != -1 and res_bc["editDistance"] < bc1_err:
                bc1_err, bc1_id = res_bc["editDistance"], bc_id_candidate

        if not bc1_id:
            SeqIO.write(record, out_files["unidentified"], "fastq")
            report.write(f"{read_id}\tAmbiguous_Barcode_5p\t{best_5p_name}\tNA\tNA\tNA\tNA\n")
            continue

        pair_3p_name, pair_3p_seq = PAIRED_3PRIME[best_5p_name]
        search_zone_3p = read_seq[best_5p_end + 1:]
        res_3p = edlib.align(pair_3p_seq, search_zone_3p, mode="HW", task="locations")
        if res_3p["editDistance"] == -1 or res_3p["editDistance"] > MAX_ERRORS_ADAPTER:
            SeqIO.write(record, out_files["unidentified"], "fastq")
            report.write(f"{read_id}\tUnidentifiable_3p\t{best_5p_name}\t{bc1_id}\tNA\tNA\tNA\n")
            continue
            
        start_3p_rel, end_3p_rel = res_3p["locations"][0]
        real_3p_start = best_5p_end + 1 + start_3p_rel
        real_3p_end = best_5p_end + 1 + end_3p_rel
        
        bc2_id, bc2_err = None, MAX_ERRORS_BARCODE + 1
        sub_read_bc2 = read_seq[max(0, real_3p_end - INDEL_BUFFER) :]
        for bc_seq, bc_id_candidate in barcode_db[pair_3p_name].items():
            res_bc = edlib.align(bc_seq, sub_read_bc2, mode="HW")
            if res_bc["editDistance"] != -1 and res_bc["editDistance"] < bc2_err:
                bc2_err, bc2_id = res_bc["editDistance"], bc_id_candidate
                
        if not bc2_id:
            SeqIO.write(record, out_files["unidentified"], "fastq")
            report.write(f"{read_id}\tAmbiguous_Barcode_3p\t{best_5p_name}\t{bc1_id}\t{pair_3p_name}\tNA\tNA\n")
            continue
            
        sample_id = sample_manifest.get(frozenset((bc1_id, bc2_id)))
        trimmed_record = record[best_5p_end + 1 : real_3p_start]

    else:
        # =====================================================================
        # БИОЛОГИЧЕСКИ КОРРЕКТНАЯ ЛОГИКА ДЛЯ SINGLE-INDEX
        # =====================================================================
        orientation = None # 'forward' или 'reverse_complement'
        ad_start, ad_end = 0, 0
        bc_id = None
        
        # 1. Сначала ищем F_raw на 5'-конце (первая половина рида)
        half_read = read_seq[:len(read_seq)//2]
        res_fwd = edlib.align(a5_seq, half_read, mode="HW", task="locations")
        
        if res_fwd["editDistance"] != -1 and res_fwd["editDistance"] <= MAX_ERRORS_ADAPTER:
            orientation = "forward"
            ad_start, ad_end = res_fwd["locations"][0]
            
            # Баркод слева от адаптера
            window_bc_end = min(len(read_seq), ad_start + INDEL_BUFFER)
            sub_read_bc = read_seq[0 : window_bc_end]
            
            bc_err = MAX_ERRORS_BARCODE + 1
            for bc_seq, candidate in barcode_db["F_raw"].items():
                res_bc = edlib.align(bc_seq, sub_read_bc, mode="HW")
                if res_bc["editDistance"] != -1 and res_bc["editDistance"] < bc_err:
                    bc_err, bc_id = res_bc["editDistance"], candidate
                    
            trimmed_record = record[ad_end + 1 :]
            
        else:
            # 2. Если на 5'-конце чисто, ищем F_revcomp на 3'-конце (вторая половина рида)
            a5_revcomp = a5_seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]
            second_half = read_seq[len(read_seq)//2:]
            res_rev = edlib.align(a5_revcomp, second_half, mode="HW", task="locations")
            
            if res_rev["editDistance"] != -1 and res_rev["editDistance"] <= MAX_ERRORS_ADAPTER:
                orientation = "reverse_complement"
                start_rel, end_rel = res_rev["locations"][0]
                # Восстанавливаем абсолютные координаты в риде
                ad_start = len(read_seq)//2 + start_rel
                ad_end = len(read_seq)//2 + end_rel
                
                # Баркод справа от адаптера (до самого конца рида)
                window_bc_start = max(0, ad_end - INDEL_BUFFER)
                sub_read_bc = read_seq[window_bc_start:]
                
                bc_err = MAX_ERRORS_BARCODE + 1
                for bc_seq, candidate in barcode_db["F_revcomp"].items():
                    res_bc = edlib.align(bc_seq, sub_read_bc, mode="HW")
                    if res_bc["editDistance"] != -1 and res_bc["editDistance"] < bc_err:
                        bc_err, bc_id = res_bc["editDistance"], candidate
                        
                trimmed_record = record[:ad_start]

        # Валидация результатов Single-indexing
        if not orientation:
            SeqIO.write(record, out_files["unidentified"], "fastq")
            report.write(f"{read_id}\tUnidentifiable_Adapter\tNA\tNA\tNA\tNA\n")
            continue
            
        if not bc_id:
            SeqIO.write(record, out_files["unidentified"], "fastq")
            report.write(f"{read_id}\tAmbiguous_Barcode\t{a5_seq}\t{orientation}\tNA\tNA\n")
            continue
            
        sample_id = sample_manifest.get(bc_id)

    # --- Общий блок записи для обоих режимов ---
    if not sample_id:
        SeqIO.write(record, out_files["unidentified"], "fastq")
        if IS_DUAL_INDEX:
            report.write(f"{read_id}\tUnknown_Combination\t{best_5p_name}\t{bc1_id}\t{pair_3p_name}\t{bc2_id}\tNA\n")
        else:
            report.write(f"{read_id}\tUnknown_Barcode\t{a5_seq}\t{orientation}\t{bc_id}\tNA\n")
        continue
        
    if sample_id not in out_files:
        out_files[sample_id] = open(os.path.join(args.outdir, f"demux_fastq/{sample_id}.fastq"), "w")
        
    SeqIO.write(trimmed_record, out_files[sample_id], "fastq")
    fully_identified += 1
    
    if IS_DUAL_INDEX:
        report.write(f"{read_id}\tFully_Identified\t{best_5p_name}\t{bc1_id}\t{pair_3p_name}\t{bc2_id}\t{sample_id}\n")
    else:
        report.write(f"{read_id}\tFully_Identified\t{a5_seq}\t{orientation}\t{bc_id}\t{sample_id}\n")

# Закрытие файлов
report.close()
for f in out_files.values():
    f.close()
print(f"Успешно завершено! Демультиплексировано ридов: {fully_identified} из {count}")