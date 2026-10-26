"""Parse the CelesTrak CSSI space-weather file (SW-All.txt) into a daily table.

Columns (FORMAT(I4,I3,I3,I5,I3,8I3,I4,8I4,I4,F4.1,I2,I4,F6.1,I2,5F6.1)):
yy mm dd BSRN ND Kp*8 Sum ap*8 Ap Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7obs Ctr81obs Lst81obs
Kp values are given x10 (e.g. 23 = 2+, 27 = 3-, 30 = 3o).
"""
import sys
import numpy as np
import pandas as pd

def parse_sw_all(path):
    rows = []
    with open(path) as f:
        in_obs = False
        for line in f:
            if line.startswith("BEGIN OBSERVED"):
                in_obs = True; continue
            if line.startswith("END OBSERVED"):
                break
            if not in_obs or not line.strip():
                continue
            # fixed-width per FORMAT(I4,I3,I3,I5,I3,8I3,I4,8I4,I4,...)
            y, m, d = int(line[0:4]), int(line[4:7]), int(line[7:10])
            kp = [int(line[18+3*i:21+3*i]) for i in range(8)]
            ksum = int(line[42:46])
            ap = [int(line[46+4*i:50+4*i]) for i in range(8)]
            Ap = int(line[78:82])
            rows.append((pd.Timestamp(y, m, d), *kp, ksum, *ap, Ap))
    cols = ["date"] + [f"Kp{i}" for i in range(8)] + ["KpSum"] + [f"ap{i}" for i in range(8)] + ["Ap"]
    df = pd.DataFrame(rows, columns=cols).set_index("date")
    return df

if __name__ == "__main__":
    df = parse_sw_all(sys.argv[1] if len(sys.argv) > 1 else "SW-All.txt")
    df.to_csv("kp_ap_daily.csv")
    print(df.index.min(), df.index.max(), len(df))
    print(df.loc["2021-10-25":"2021-11-01"])
    # sanity: Ap should equal round(mean(ap0..ap7))
    chk = (df[[f"ap{i}" for i in range(8)]].mean(axis=1).round() - df["Ap"]).abs()
    print("max |mean(ap)-Ap| =", chk.max(), " rows>1:", int((chk > 1).sum()))
