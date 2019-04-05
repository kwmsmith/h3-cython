from skbuild import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="h3cy",
    version="0.0.3",
    description="Python bindings (via Cython) to the H3 hierarchical geospatial indexing library.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kwmsmith/h3-cython",
    author="Kurt Smith",
    email="kwmsmith@gmail.com",
    packages=["h3cy"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: C",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
    ],
)
