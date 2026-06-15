import os
import subprocess
import msprime
import pyslim
import numpy as np

POP_SIZE = 10_000 #ancestral N (N0)
SEQ_LEN = 1_000_000
BASE_MU = 2.5e-7        
RECOMB_RATE = 1e-8      
N_REPS = 50             
TAU_MAX = 8.0 #target scaled time (generations/2N)

#decline ratios 
RATIOS = [1.0/4.0, 1.0/16.0, 1.0/64.0]

OUT_DIR = "results_decline_n0_10000"
os.makedirs(OUT_DIR, exist_ok=True)

#burn-in phase (msprime)
def run_burnin(rep):
    seed = 21 + rep
    out_trees = f"{OUT_DIR}/burnin_rep{rep}.trees"
    
    if os.path.exists(out_trees):
        return out_trees

    print(f"\n=== Burn-in replicate {rep} ===")
    demography = msprime.Demography()
    demography.add_population(initial_size=POP_SIZE)

    ts = msprime.sim_ancestry(
        samples=POP_SIZE,
        ploidy=2,
        demography=demography,
        random_seed=seed,
        recombination_rate=RECOMB_RATE,
        sequence_length=SEQ_LEN
    )

    ats = pyslim.annotate(ts, model_type="WF", tick=1, stage="late")
    
    #add mutations for standing variation
    next_id = pyslim.next_slim_mutation_id(ats)
    mut_model = msprime.SLiMMutationModel(type=1, next_id=next_id)
    vts = msprime.sim_mutations(
        ats,
        rate=BASE_MU,
        model=mut_model,
        keep=True,
        random_seed=seed + 1
    )

    vts.dump(out_trees)
    return out_trees

#run SLiM
def run_slim(rep, ratio, burnin_file):
    end_n = int(POP_SIZE * ratio)
    num_gen = int(TAU_MAX * 2 * end_n)
    sample_interval = int(max(10, num_gen / 200))
    current_mu = BASE_MU
    
    #output filename
    ratio_name = f"{ratio:.4f}".replace(".", "p")
    out_csv = f"{OUT_DIR}/pi_ratio{ratio_name}.csv"
    
    #create SLiM Command
    cmd = [
        "slim",
        "-d", f'INIT_TREES=\\"{burnin_file}\\"',
        "-d", f"MU={current_mu}",
        "-d", f"START_N={POP_SIZE}",
        "-d", f"END_N={end_n}",
        "-d", f"NUM_GEN={num_gen}",
        "-d", f"SAMPLE_INTERVAL={sample_interval}",
        "-d", f"REP={rep}",
        "-d", f'OUTFILE=\\"{out_csv}\\"',
        "decline_instant.slim"
    ]

    subprocess.run(" ".join(cmd), shell=True, check=True)

#main workflow
if __name__ == "__main__":
    for rep in range(1, N_REPS + 1):
        burnin_file = run_burnin(rep)
        for ratio in RATIOS:
            run_slim(rep, ratio, burnin_file=burnin_file)

    print("\nAll simulations complete")