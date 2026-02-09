import sqlite3
import os

db_path = "pksmart.db"

if not os.path.exists(db_path):
    print(f"Database {db_path} not found. Skipping migration.")
    exit(0)


def add_column():
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    tables = ["projects", "ind_reports"]
    for table in tables:
        try:
            cursor.execute(
                f"ALTER TABLE {table} ADD COLUMN is_deleted BOOLEAN DEFAULT 0"
            )
            print(f"Added is_deleted column to {table}")
        except sqlite3.OperationalError as e:
            if "duplicate column name" in str(e):
                print(f"Column is_deleted already exists in {table}")
            else:
                print(f"Error altering {table}: {e}")

    conn.commit()
    conn.close()


if __name__ == "__main__":
    add_column()
