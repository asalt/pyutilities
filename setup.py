from setuptools import setup
def calculate_version(inputfile):
    version_list =  [x.split('\'')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

package_version = calculate_version('pyutilities.py')

setup(
    name='PyUtilities',
    version=package_version,
    py_modules=['pyutilities'],
    install_requires=[
        'Click', 'beautifulsoup4', 'lxml', 'biopython',
    ],
    entry_points="""
    [console_scripts]
    pyutilities=pyutilities:cli
    """
)
