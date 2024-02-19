import subprocess

def check_version(command):
    try:
        result = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT).decode('utf-8').strip()
        return result
    except subprocess.CalledProcessError as e:
        return e.output.decode('utf-8').strip()

def main():
    software_commands = {
        'Python': 'python --version',
        'Conda': '/opt/conda/bin/conda --version',
        'Mamba': '/opt/conda/bin/mamba --version',
        'Busco': 'busco --version',
        'Muscle': 'muscle -version',
        'MAFFT': 'mafft --version',
        'Trimal': 'trimal --version',
        'IQ-TREE': 'iqtree --version'
    }

    for program, command in software_commands.items():
        version_message = check_version(command)
        print(f"{program}: {version_message}")

if __name__ == "__main__":
    main()