# Changelog
All notable changes to this project will be documented in this file. See [conventional commits](https://www.conventionalcommits.org/) for commit guidelines.

- - -
## [v0.4.0](https://github.com/yoctoyotta1024/SDM-1D-KiD/compare/2baae7441fbe53a16fa13447190d6bd3a0284470..v0.4.0) - 2026-11-25
#### Bug Fixes
- missing logger arg - ([28de19d](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/28de19d60c58e54ce0847535a057736c75eab12e)) - clara.bayley
- wrong name for bin edges file - ([f80761c](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/f80761cf38f36abd77adb75cbde6c7b20b360f62)) - clara.bayley
- correct axes - ([804c8cd](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/804c8cd68a89d22e696afe5f392ee58edda8541a)) - clara.bayley
- path missing - ([24e1ce4](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/24e1ce45aa8e5b779ad4d41b3f32633e767a16e3)) - clara.bayley
#### Continuous Integration
- don't run model twice - ([214776d](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/214776dc452e84cd46ffb1764613c6dcad05c7fc)) - clara.bayley
#### Documentation
- include post-processing pysdm information - ([e45133d](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/e45133d3cf7ba60d87c51edda214ec03bbfb8feb)) - clara.bayley
- include post-processing pysdm information - ([635a984](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/635a98470318c125ae081670113ea0946b802f3d)) - clara.bayley
- document on cleo vs pysdm differences - ([f6fb1e8](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/f6fb1e878b1974b44c78ed98816e7d950e927676)) - clara.bayley
#### Features
- implementation for post-processing PySDM output - ([3f3e4cf](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/3f3e4cf76750b44e42e21d72b13bd1fd511275bc)) - clara.bayley
- use alpha sampling in PySDM - ([e300e22](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/e300e22339068ccf2c3e644bfd28115361a4d13e)) - clara.bayley
- implementation for running PySDM - ([5c6d1fc](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/5c6d1fc24623e08ebfa910b1d09098c5b2403fed)) - clara.bayley
- implementation for inputfiles for PySDM - ([d680acd](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/d680acd714d516063e2bccf051225979b37c2966)) - clara.bayley
- placeholder files for PySDM 1-D KiD - ([2baae74](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/2baae7441fbe53a16fa13447190d6bd3a0284470)) - clara.bayley
#### Miscellaneous Chores
- rename CLEO Cleo - ([d2370a0](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/d2370a0894491bcb766eff10a567659f15ce2517)) - clara.bayley
- rename CLEO Cleo - ([f9108ac](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/f9108ac193329cae3ab3a09ef7b4967ac2e5fcff)) - clara.bayley
- move common.py outside pympdata_driver - ([cae3b2c](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/cae3b2c0f627043be5d0d76df9a5a7e662ead3c0)) - clara.bayley
#### Performance Improvements
- correct docstring - ([eeba327](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/eeba3271495924a67b127141a236076542ab89c5)) - clara.bayley
- bulk scheme only needs 1 ensemble member - ([afea5d5](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/afea5d5e1f26d239f5a502317087bc74bad4521b)) - clara.bayley
- format progress counter - ([c591cd7](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/c591cd7dbf73041bada2a907a43333b6a221d675)) - clara.bayley
- log when call completes - ([4d6a49e](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/4d6a49e05cc3b221ef59ade35d01b1c9a21fea93)) - clara.bayley
- minor renaming - ([7a7ffe6](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/7a7ffe6b02a6ede384adc00ba8706af30805ec65)) - clara.bayley
- correct date - ([d2ea75b](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/d2ea75b04757ed2d5a722e39777c44f4e3967ac7)) - clara.bayley
- safe read floats from config - ([ce033e8](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/ce033e863cf60ce85efd0785a0caa2b2ec0554c0)) - clara.bayley
- corrections - ([0b019ae](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/0b019aef1ff66c8facccca5ac93f9766fb0bf741)) - clara.bayley
#### Refactoring
- add post-processing for products to make variables in Hill plots - ([4de9f97](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/4de9f977c1e5f6992002c9dc9662e4c8fd196658)) - clara.bayley
- move add_radar_reflectivity into common - ([90afcde](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/90afcde8dee949d9f1a77742563d1d2960641dec)) - clara.bayley
- add products to PySDM output required for Hill plots - ([e82ac7a](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/e82ac7a210f8fee4557b49ef6f6d098c7aa5f087)) - clara.bayley
- match cleo and pysdm metadata better - ([164b23a](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/164b23adc202674e87de22ffdd73625216c094c4)) - clara.bayley
- use kappa=1.243 for PySDM - ([8cdfd95](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/8cdfd9539b7f129d5a96004c7bee14530354b185)) - clara.bayley
- add msol distribution to pysdm post-processing - ([75e3a80](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/75e3a806401a08c3e3825618d8f2f0fe71d8f579)) - clara.bayley
- post-process msol to number and mass of dry aersol distrib in Cleo - ([63ae49c](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/63ae49c9d3bb4ec6165ac161af7e7a8ca0506f0f)) - clara.bayley
- move common post-processing functions into common - ([16befd7](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/16befd71f1a1b694a2f3fdeee1c1fb4816e79bc1)) - clara.bayley
- minor changes to Cleo post-processing - ([5ca038d](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/5ca038d6c0142750001b10f61601c51f654b9ca3)) - clara.bayley
- minor changes to match PyMPDATA to PySDM - ([5501917](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/550191747494d6b650501c482faffd2780009f85)) - clara.bayley
- allow PySDM to use Exner thermodynamic profiles - ([d46dee4](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/d46dee46341b6aceef749d2779aca5aa4af9daba)) - clara.bayley
- allow PySDM to take Hydrodynamic or Long collision kernel - ([060f388](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/060f38880947116abaa684b4915804437386ac5b)) - clara.bayley
- make PySDM KiD setup similar or identical to microphysics schemes run directly with PYMPDATA - ([133b37a](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/133b37a831e3afb3fb3c85f54b704703c85db58e)) - clara.bayley
- add print statement to pysdm runtime - ([4f65c82](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/4f65c82a9e2269d6efad1ee8b310cb263955c75e)) - clara.bayley

- - -

## [v0.3.1](https://github.com/yoctoyotta1024/SDM-1D-KiD/compare/924ebe507f9dacff83d3613d01fb7d8ba497d695..v0.3.1) - 2026-11-20
#### Bug Fixes
- use random seeds for alpha sampling - ([a6f0513](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/a6f05130773ed9fd78779eeeec286aa4626fcab3)) - clara.bayley
#### Miscellaneous Chores
- use alpha=1.0 and xi_by_pressure=True - ([fc8b49b](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/fc8b49b5adb893b274e388920428414df47d6d0f)) - clara.bayley
#### Refactoring
- enable xi_by_pressure initialisation - ([924ebe5](https://github.com/yoctoyotta1024/SDM-1D-KiD/commit/924ebe507f9dacff83d3613d01fb7d8ba497d695)) - clara.bayley

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
