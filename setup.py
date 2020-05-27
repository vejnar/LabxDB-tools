import shutil
from distutils.core import setup
from distutils.command.install_scripts import install_scripts

class remove_extension_install_scripts(install_scripts):
    def run(self):
        install_scripts.run(self)
        for script in self.get_outputs():
            if script.endswith('.py'):
                shutil.move(script, script[:-3])

setup(
    name='labxdb',
    version='1.0',
    packages=['labxdb'],
    author='Charles E. Vejnar',
    url='https://gitlab.com/vejnar/labxdb-python',
    license='Mozilla Public License 2.0 (MPL 2.0)',
    long_description=open('README.md').read(),
    install_requires=['pyfnutils', 'requests', 'urllib3'],
    cmdclass = {'install_scripts': remove_extension_install_scripts},
    scripts=['script/export_sra_fastq.py', 'script/export_sra.py', 'script/fastq_strip.py', 'script/import_seq.py', 'script/import_sra_id.py', 'script/import_sra.py']
)
