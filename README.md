#STRAssembly

This is a software tool relying on GATB-CORE library.

The architecture of the tool is as follows:

    * a CMakeLists.txt file used for building the project
    * a 'tools' directory holding a default source code using GATB-Core
    * a 'scripts' directory holding a script to automatically package the tool
    * a 'thirdparty' directory holding the gatb-core resources
    * a 'doc' directory
    * a 'tests' directory holding test procedures

The 'thirdparty' directory is only available for tool created outside the GATB-Tools repository.
Tools located within GATB-Tools rely on a common GATB-Core sub-module already available in this repository.

It is advised to use:

    * 'tests' directory to hold test procedures: scripts and/or small sized data files
    * 'scripts' directory to hold any scripts this tool relies on
    * 'doc' directory to hold tool's documentation
    * 'thirdparty' directory to hold any third-party librairies this tool relies on

It is worth noting that 'tools' directory is organised using sub-folders; by default, there is
at least one such sub-folder called 'STRAssembly'. It holds the source code of the tool. However, when
considering a more complex software, it could be nice to setup several "sub-tools", each of them
making a particular feature. In that case, you can easily create several "sub-tool" folders inside
"tools", each of them having a "src" folder containing the source code, as well as a "main.cpp", for
each feature. Using this organisation has a big advantage: the provided CMakeLists.txt is aware of
that, so you do not have to edit the CMake file when you add a new "sub-tool". As a real example, you
can have a look at the DiscoSNP software.

#License

Please not that GATB-Core is distributed under Affero-GPL license.

#Dependencies

The following third parties should be already installed:

* cmake 2.8+ (mandatory)

#Download

Clone the project with gatb-core.

```
git clone --recursive git@github.com:gulfemd/STRAssembly
```

#Project build

For building your project, you should do the following

```
mkdir build
cd build
cmake ..
make -j8
```

Then, you should get a binary holding the name of the project within 'build/tools'.

Note: the first compilation should take some time since the GATB-CORE library is generated.

#Usage

```
./build/tools/STRAssembly -out <OUTPUT_DIR> -sam <READS_FILE_PATH> -nb-cores 0 -verbose 1
```
