from setuptools import setup, find_packages

setup(
    name="GEx_Predictor",
    version="1.0.0",
    author="AdriaVico",
    author_email="adria.vico@irbbarcelona.org",
    description="Predict gene expression from SMILES using GLOBAL molecular signatures.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Adriv020/GEx_Predictor",
    packages=find_packages(),
    install_requires=[
        "torch>=2.0.0,<3.0.0",
        "rdkit-pypi>=2023.3.1b1",
        "signaturizer>=0.1.0",
        "certifi>=2024.0.0",
        "numpy>=1.24.0",
    ],
    python_requires=">=3.10",
    license="MIT",
    keywords="gene expression SMILES chemistry deep learning bioinformatics",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)