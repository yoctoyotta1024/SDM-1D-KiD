# Changelog
All notable changes to this project will be documented in this file. See [conventional commits](https://www.conventionalcommits.org/) for commit guidelines.

- - -
## [v0.3.0](https://github.com/yoctoyotta1024/SDM-1D-KiD/compare/2c42884302b356f069305844c2c67fb9bb9b18ec..v0.3.0) - 2026-11-20
#### Bug Fixes
- drop rhod if similar to rho_dry for bulk scheme - ([f29c6a9](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/f29c6a9373e6c201db103607e5f14c12b416cc19)) - clara.bayley
- use json.loads not dict - ([cbd2d88](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/cbd2d88941e3eaf5a195b6ffcc2a215573083376)) - clara.bayley
- launch cleo runs from seperate MPI processes - ([8ef2a2e](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/8ef2a2e1bbad7f83965d83e3ff85075c4ccc946d)) - clara.bayley
- incorrect import - ([c4156dc](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/c4156dca58fa446ff78b1e409bd589848a193db9)) - clara.bayley
- missing allsupers from send_dynamics - ([030b478](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/030b4781267ed190822b9d4b38086038f9f72676)) - clara.bayley
#### Continuous Integration
- add bulk post-processing to ci and change ensemble labels - ([48243df](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/48243dfbf874df58c356a9246da5fc912e7551c8)) - clara.bayley
- add running of Cleo 1-D KiD to CI - ([a90a974](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/a90a974665b16917ff53dc2ce135ac267fddcd07)) - clara.bayley
#### Documentation
- minor corrections - ([2741c59](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/2741c59246d57434c4d1b4b099aa57b860e4b91c)) - clara.bayley
- add details of bulk scheme post-processing - ([9d0ff87](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/9d0ff8761447af04a98688f46728c9de1944ca89)) - clara.bayley
- docs on post-processing cleo output - ([a91ab58](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/a91ab58a21e57a1069de0eccaef1ee8a65a9ed0a)) - clara.bayley
- fix multiple references - ([be97a34](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/be97a345bbdd9b4ba1129585f73632ee256a426f)) - clara.bayley
- more info on Cleo runs - ([c546fe7](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/c546fe7e1f5d7d377a78fd42579c2155b1c89544)) - clara.bayley
- minor re-formatting - ([2c42884](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/2c42884302b356f069305844c2c67fb9bb9b18ec)) - clara.bayley
#### Features
- add script for postprocessing bulk scheme - ([560ac47](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/560ac47c3c8894fe610b9ae2a9b9c9b8d6f31b93)) - clara.bayley
- add postprocesssing of superdroplets for Hill paper variables - ([37ef38b](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/37ef38b6fbd5f51e80da04748a9a6cb12c95c874)) - clara.bayley
- add script for postprocessing Cleo output - ([29e21ba](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/29e21ba8d88e5b041a9122c082fe4116098b601a)) - clara.bayley
- scripts to run PyMPDATA with bulk scheme - ([b02c687](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/b02c687b9c5cd2331dfc9b864f822f527ba85f00)) - clara.bayley
- files to run 1-D KiD using PyMPDATA - ([ef497ca](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/ef497cabac819887679abbb6d510197154d0fe8c)) - clara.bayley
- bash script to call run_cleo.py - ([5379b9a](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/5379b9a595b884e09bf4955392afec66377376a5)) - clara.bayley
- files to run Cleo coupled to arbitrary numpy thermodynamics class - ([7db899d](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/7db899d3a9e89162f35c7c8bbb239fba5b6e5cad)) - clara.bayley
- files for Cleo initial conditions - ([ff3b591](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/ff3b5916a69950d851d0dc03acef8c22373a4271)) - clara.bayley
- files to compile Cleo - ([99b72a0](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/99b72a010edf4c061c905fcb02ef501a43f4ddb5)) - clara.bayley
- placeholder files for Cleo 1-D KiD - ([c7d3368](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/c7d336839b9273eb6be9d8f96e23911871f4cb1e)) - clara.bayley
#### Miscellaneous Chores
- rename post_process to postprocess - ([d38b37e](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/d38b37e1f2bfa1e6aa8fd5e27a5a13ca6ee3ffac)) - clara.bayley
- increase memory allowance to run Cleo through SLURM - ([89c87ea](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/89c87ea07d6f0af90f8794d7fac801d90527135f)) - clara.bayley
- spelling - ([a45244b](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/a45244bfdb23deece5c41efe29adcad6a424b693)) - clara.bayley
- move thermo module into pympdata_driver module - ([83c878f](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/83c878ffc920efdeb1c6fd82800eaee60950c3bf)) - clara.bayley
- correct created date - ([3f78b18](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/3f78b1872771dd299974a817d0b311e5ee0a7573)) - clara.bayley
#### Performance Improvements
- faster loading when python path included needed directories - ([b34ee60](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/b34ee60a5428bcf8f9d05d52baee67a2fa75f63d)) - clara.bayley
- delete if file exisits before over-writing - ([e5167ee](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/e5167ee50585d25575e6938ef9c7970b2bee192c)) - clara.bayley
- capitalise print statement - ([8e65e37](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/8e65e370d559fbf4dd0ff4481c8c7283b9720c08)) - clara.bayley
#### Refactoring
- set xi_by_pressure in cleo config - ([49de838](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/49de83896d8a54173e9c5c598705c5437d58d379)) - clara.bayley
- move KiD parameters into run_cleo and run_bulk_scheme - ([8f89223](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/8f89223b74ae9968d5d0d8e3eb7bd8885a242baf)) - clara.bayley
- match bulk scheme metadata of output to cleo postprocessing - ([603f78a](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/603f78a90c3f744b93bfd8b7e0bbc2be43efd624)) - clara.bayley
- be clearerin cleo metadata about solute included/excluded in calculations - ([4407e78](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/4407e787257b20dd88d0d50a420efb439fc3e94b)) - clara.bayley
- add more metadata to cleo output - ([693e986](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/693e9869e1a67c16a124ae66e98bc314df942d67)) - clara.bayley
- move geomean and geosig into run_cleo parameters - ([ebad3d7](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/ebad3d79ebf639ab6f6bac83b6e8a53f267a4830)) - clara.bayley
- allow bulk scheme runs to output data - ([4c8f5e6](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/4c8f5e624ed79d1289a37b97c74637a54f95acab)) - clara.bayley
- add SLURM options to run SDM bash - ([d414169](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/d414169db2643e61f199d9a0512ec427126eb2e1)) - clara.bayley
- enable ensemble of runs creation with run_cleo - ([a81227d](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/a81227dcf43c1116f397e12d8745259d3a2fe3d4)) - clara.bayley

- - -

## [v0.2.0](https://github.com/yoctoyotta1024/SDM-1D-KiD/compare/8bf362d7915b41717c034e858f6287058e9341b6..v0.2.0) - 2026-11-17
#### Documentation
- update get-started and contributing guidelines - ([553b803](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/553b803290ce3aae5069b2a6038964df3640ef35)) - clara.bayley
#### Features
- add clang-format to base repo - ([8bf362d](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/8bf362d7915b41717c034e858f6287058e9341b6)) - clara.bayley
#### Performance Improvements
- missing space in readme - ([fec2f24](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/fec2f244cdefdac1f3bf286039346376a8612231)) - clara.bayley

- - -

## [v0.1.0](https://github.com/yoctoyotta1024/SDM-1D-KiD/compare/23821afcbd4ac0755e7df2e2b9681861083eaaf8..v0.1.0) - 2026-11-17
#### Features
- add basic repository boiler plate - ([fec1448](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/fec1448b75dceb6a1d3a9cee2b0d2126c3a81fe3)) - clara.bayley
- add basic documentation - ([23821af](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/23821afcbd4ac0755e7df2e2b9681861083eaaf8)) - clara.bayley

- - -

Changelog generated by [cocogitto](https://github.com/cocogitto/cocogitto).
