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
        "bin/clear_data.py",
        "bin/convertPIMCuuid.py",
        "bin/convert_state.py",
        "bin/energyfit.py",
        "bin/fix_tiny.py",
        "bin/gensubmit.py",
        "bin/getresults.py",
        "bin/hescrew.py",
        "bin/make_pbs_resubmit_script.py",
        "bin/make_pbs_script.py",
        "bin/merge.py",
        "bin/mergeBK.py",
        "bin/numform.py",
        "bin/pimcbin.py",
        "bin/plot_matrix_estimator.py",
        "bin/plot_scalar_estimator.py",
        "bin/plot_tensor_estimator.py",
        "bin/plot_vector_estimator.py",
        "bin/plotbin.py",
        "bin/plotoptions.py",
        "bin/procisf.py",
        "bin/rename.py",
        "bin/rmempty.py",
        "bin/rsubmit.py",
        "bin/timeStepScaling.py",
        "bin/wikiUpload.py",
        ],
)
