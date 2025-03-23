import mysql.connector
from mysql.connector import Error

def add_protein_sequence_column():
    """Adds the 'protein_sequence' column to the existing table if it doesn't exist."""
    try:
        conn = mysql.connector.connect(**MYSQL_CONFIG)
        cursor = conn.cursor()

        # Check if the 'protein_sequence' column exists in the 'bioactivity' table
        cursor.execute("""
        SELECT COUNT(*) 
        FROM INFORMATION_SCHEMA.COLUMNS 
        WHERE table_name = 'bioactivity' AND column_name = 'protein_sequence';
        """)

        # If the column doesn't exist, add it
        if cursor.fetchone()[0] == 0:
            cursor.execute("""
            ALTER TABLE bioactivity 
            ADD COLUMN protein_sequence TEXT;
            """)
            conn.commit()
            print("✅ Column 'protein_sequence' added to the table!")
        else:
            print("✅ Column 'protein_sequence' already exists!")

        cursor.close()
        conn.close()
    except Error as e:
        print(f"❌ Error adding column to table: {e}")