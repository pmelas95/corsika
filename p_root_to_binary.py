import uproot
import numpy as np
import sqlite3
import argparse

def extract_from_root(root_file, tree_name="gRooTracker2"):
    with uproot.open(root_file) as file:
        tree = file[tree_name]

        # Extract branches
        pdg = np.concatenate(tree["StdHepPdg"].array(library="np"))
        n_particles = tree["StdHepN"].array(library="np")
        p4 = np.concatenate(tree["StdHepP4"].array(library="np"))
        pos = np.concatenate(tree["EvtVtx"].array(library="np"))

        # New metadata fields
        has_min_dt = np.concatenate(tree["Has_MIN_DT"].array(library="np"))
        has_min = np.concatenate(tree["Has_MIN"].array(library="np"))
        has_dt = np.concatenate(tree["Has_DT"].array(library="np"))
        has_dt_min = np.concatenate(tree["Has_DT_MIN"].array(library="np"))
        has_y = np.concatenate(tree["Has_Y"].array(library="np"))
        has_gene = np.concatenate(tree["Has_Gene"].array(library="np"))

        # Initialize flat arrays
        px, py, pz, e = [], [], [], []
        x, z, t, shower_ids = [], [], [], []

        particle_index = 0
        for shower_id, n in enumerate(n_particles):
            x.append(pos[shower_id * 4 + 0])
            t.append(pos[shower_id * 4 + 1])
            z.append(pos[shower_id * 4 + 2])

            for i in range(n):
                px.append(p4[particle_index][0])
                py.append(p4[particle_index][1])
                pz.append(p4[particle_index][2])
                e.append(p4[particle_index][3])
                shower_ids.append(shower_id)
                particle_index += 1

        return (shower_ids, pdg, px, py, pz, x, z, t, e,
                n_particles, has_min_dt, has_min, has_dt, has_dt_min, has_y, has_gene)

def save_to_sqlite(db_file, shower_ids, pdg, px, py, pz, x, z, t, e,
                   has_min_dt, has_min, has_dt, has_dt_min, has_y, has_gene):
    with sqlite3.connect(db_file, timeout=10) as conn:
        cursor = conn.cursor()

        # Create tables
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS particles (
                shower INT,
                pdg INT,
                px FLOAT,
                py FLOAT,
                pz FLOAT,
                x FLOAT,
                z FLOAT,
                t FLOAT,
                e FLOAT,
                has_min_dt FLOAT,
                has_min FLOAT,
                has_dt FLOAT,
                has_dt_min FLOAT,
                has_y FLOAT,
                has_gene FLOAT
            )
        ''')

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS showers (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                nparticles INTEGER
            )
        ''')

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS input (
                nparticles INTEGER,
                version REAL,
                nshow INTEGER,
                nshowsim INTEGER,
                model_high INTEGER,
                model_low INTEGER,
                eslope REAL,
                erange_high REAL,
                erange_low REAL,
                ecuts_hadron REAL,
                ecuts_muon REAL,
                ecuts_electron REAL,
                ecuts_photon REAL
            )
        ''')

        # Insert particle records
        for i in range(len(pdg)):
            shower_id = shower_ids[i]
            if pdg[i] == 0:
                continue

            cursor.execute('''
                INSERT INTO particles (
                    shower, pdg, px, py, pz, x, z, t, e,
                    has_min_dt, has_min, has_dt, has_dt_min, has_y, has_gene
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                shower_id,
                int(pdg[i]),
                float(px[i]), float(py[i]), float(pz[i]),
                float(x[shower_id]), float(z[shower_id]), float(t[shower_id]), float(e[i]),
                float(has_min_dt[i]),
                float(has_min[i]),
                float(has_dt[i]),
                float(has_dt_min[i]),
                float(has_y[i]),
                float(has_gene[i])
            ))

        conn.commit()

        # Shower table
        insert_shower_data(conn)

        # Input config table
        cursor.execute("SELECT COUNT(*) FROM particles")
        nparticles = cursor.fetchone()[0]
        parse_and_insert_input(conn, nparticles)

# Function to insert shower data into the database
def insert_shower_data(conn):
    start_id = 0
    max_id = 1232129

    with conn:
        cursor = conn.cursor()

        for current_id in range(start_id, max_id + 1):
            cursor.execute('INSERT INTO showers (id, nparticles) VALUES (?, ?)', (current_id, 0))

        conn.commit()

# Function to insert input data
def parse_and_insert_input(conn, nparticles):
    version = 7.7
    run_number = 1
    number_of_showers = 1232130
    model_high = 3
    model_low = 3
    energy_slope = -2.7
    energy_min = 1.3
    energy_max = 100000 
    ecuts = [0.05, 0.05, 0.05, 0.05]

    with conn:
        cursor = conn.cursor()
        cursor.execute('''
            INSERT INTO input (
                nparticles, version, nshow, nshowsim,
                model_high, model_low, eslope,
                erange_high, erange_low,
                ecuts_hadron, ecuts_muon, ecuts_electron, ecuts_photon
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            nparticles, version, number_of_showers, number_of_showers,
            model_high, model_low, energy_slope,
            energy_max, energy_min, *ecuts
        ))
        conn.commit()
    print("Inserted input data.")

def root_to_db(root_file, db_file):
    (shower_ids, pdg, px, py, pz, x, z, t, e,
     n_particles, has_min_dt, has_min, has_dt, has_dt_min, has_y, has_gene) = extract_from_root(root_file)

    save_to_sqlite(db_file, shower_ids, pdg, px, py, pz, x, z, t, e,
                   has_min_dt, has_min, has_dt, has_dt_min, has_y, has_gene)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert ROOT file to SQLite with muon bundle metadata.")
    parser.add_argument('--root_file', required=True, help="Path to the gRooTracker2 ROOT file")
    parser.add_argument('--db_name', default='corsika_simulation.db', help="Output SQLite DB name")
    args = parser.parse_args()

    root_to_db(args.root_file, args.db_name)
