import numpy as np
import matplotlib.pyplot as plt

# ─── Paramètres ───────────────────────────────────────────────────────────────
TRAJECTORY_FILE = "trajectory.xyz"   # chemin vers ton fichier
OUTPUT_PNG      = "msd.png"
# ──────────────────────────────────────────────────────────────────────────────


def parse_trajectory(filename):
    """Lit le fichier .xyz et retourne positions, ids molécules, boîte."""
    frames = []
    mol_ids_all = []
    box_sizes = []

    with open(filename) as f:
        while True:
            # Ligne 1 : nombre d'atomes
            line = f.readline()
            if not line:
                break
            n_atoms = int(line.strip())

            # Ligne 2 : header (cell size, step, ...)
            header = f.readline().strip()
            # Extraire la taille de la boîte
            cell_str = [tok for tok in header.split() if tok.startswith("cell:")]
            if cell_str:
                cell_vals = cell_str[0].split(":")[1].split(",")
                box = np.array([float(v) for v in cell_vals])
            else:
                box = None

            # Lire les atomes
            positions = []
            mol_ids = []
            for _ in range(n_atoms):
                parts = f.readline().split()
                mol_ids.append(int(parts[0]))
                # parts[1] = species, parts[2:5] = x y z
                positions.append([float(parts[2]), float(parts[3]), float(parts[4])])

            frames.append(np.array(positions))       # (N, 3)
            mol_ids_all.append(np.array(mol_ids))    # (N,)
            box_sizes.append(box)                    # (3,)

    return frames, mol_ids_all, box_sizes


def unwrap_trajectory(frames, box_sizes):
    """Déplie les trajectoires pour corriger les sauts PBC."""
    n_frames = len(frames)
    n_atoms  = frames[0].shape[0]
    unwrapped = np.zeros((n_frames, n_atoms, 3))
    unwrapped[0] = frames[0]

    for t in range(1, n_frames):
        box = box_sizes[t]
        delta = frames[t] - frames[t - 1]
        # Correction PBC : ramène le déplacement dans [-box/2, box/2]
        delta -= box * np.round(delta / box)
        unwrapped[t] = unwrapped[t - 1] + delta

    return unwrapped  # (n_frames, n_atoms, 3)


def compute_msd_atoms(unwrapped):
    """MSD moyen sur toutes les particules."""
    r0 = unwrapped[0]  # positions initiales
    msd = np.mean(np.sum((unwrapped - r0) ** 2, axis=2), axis=1)
    return msd


def compute_msd_com(unwrapped, mol_ids):
    """MSD des centres de masse des molécules."""
    unique_mols = np.unique(mol_ids)
    n_frames = unwrapped.shape[0]
    com_traj = np.zeros((n_frames, len(unique_mols), 3))

    for i, mol in enumerate(unique_mols):
        mask = mol_ids == mol
        com_traj[:, i, :] = unwrapped[:, mask, :].mean(axis=1)

    r0 = com_traj[0]
    msd = np.mean(np.sum((com_traj - r0) ** 2, axis=2), axis=1)
    return msd


def main():
    print("Lecture de la trajectoire...")
    frames, mol_ids_all, box_sizes = parse_trajectory(TRAJECTORY_FILE)
    n_frames = len(frames)
    print(f"  {n_frames} frames, {frames[0].shape[0]} atomes par frame")

    # On utilise les mol_ids du premier frame (constants)
    mol_ids = mol_ids_all[0]

    print("Dépliage PBC...")
    unwrapped = unwrap_trajectory(frames, box_sizes)

    print("Calcul MSD particules...")
    msd_atoms = compute_msd_atoms(unwrapped)

    print("Calcul MSD centres de masse...")
    msd_com = compute_msd_com(unwrapped, mol_ids)

    # Axe temporel (en steps)
    steps = np.arange(n_frames)

    # ─── Plot ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(steps, msd_atoms, label="Particules (moyenne)", linewidth=1.5)
    ax.plot(steps, msd_com,   label="Centres de masse mol.", linewidth=1.5, linestyle="--")
    ax.set_xlabel("Pas de temps")
    ax.set_ylabel("MSD")
    ax.set_title("Mean Square Displacement")
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=150)
    print(f"Plot sauvegardé : {OUTPUT_PNG}")
    plt.show()

    # ─── Sauvegarde des données ────────────────────────────────────────────
    np.savetxt("msd_atoms.dat", np.column_stack([steps, msd_atoms]),
               header="step  MSD_atoms")
    np.savetxt("msd_com.dat",   np.column_stack([steps, msd_com]),
               header="step  MSD_com")
    print("Données sauvegardées : msd_atoms.dat, msd_com.dat")


if __name__ == "__main__":
    main()
