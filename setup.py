import pathlib
import setuptools

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setuptools.setup(
    name="pimcscripts",
    version="0.1",
    packages=setuptools.find_packages(),
    description="Scripts for analzying the results of path integral quantum Monte Carlo simulations",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/DelMaestroGroup/pimcscripts",
    author="Adrian Del Maestro",
    author_email="adrian@delmaestro.org",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    include_package_data=True,
    install_requires=["numpy", "matplotlib", "scipy", "docopt", "joblib"],
    scripts=["bin/pimcave.py",
        "bin/pimcplot.py",
        "bin/reduce-one.py",
        ],
)

# "gt_cotccp=graphenetools.gt_c_one_third_commensurate_command_plot:main",
