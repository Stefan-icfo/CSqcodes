from qcodes import initialise_or_create_database_at
import subprocess
import os


#DATABASE_LOCATION = ".\Data\Raw_data\QuantumSimulator.db"
#DATABASE_LOCATION = ".\Data\Raw_data\CD11_D7_C1_zurichdata.db"
DATABASE_LOCATION = ".\Data\Raw_data\CD13_E3_C2.db"
#DATABASE_LOCATION = ".\Data\Raw_data\CD20_f2top.db"
#DATABASE_LOCATION = ".\Data\Raw_data\testruns_withdevice.db"
#DATABASE_LOCATION = ".\Data\Raw_data\testruns_nodevice.db"
BACKUP_DATABASE_LOCATION = "Z:"+"\\"+"Electromechanics"+"\\"+"Projects"+"\\"+"chargesensor"+"\\"+"backups"+"\\"+"Raw_data"+"\\"+'CD13_E3_C2_backup.db'

#TEST_DATABASE_LOCATION=

initialise_or_create_database_at(DATABASE_LOCATION)
#initialise_or_create_database_at("C:\Users\sforstner\Desktop\testdb.db")
#initialise_or_create_database_at(BACKUP_DATABASE_LOCATION)
#BACKUP_NAME = f'CD11_D7_C1_backup.db'

#initialise_or_create_database_at(DATABASE_LOCATION)


#def backupDatabase():
#    subprocess.run(['sqlite3', f"{DATABASE_LOCATION}",
#                    f".backup '{BACKUP_DATABASE_LOCATION}'"])

def backupDatabase():
    command = f".backup '{BACKUP_DATABASE_LOCATION}'"
    subprocess.run(['sqlite3', DATABASE_LOCATION], input=command, text=True)