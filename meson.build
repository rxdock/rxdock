project('RxDock', 'cpp', version : '0.1.0', default_options : ['cpp_std=c++11'])

incRbt = include_directories('include')
install_headers(
  files('include/rxdock/AlignTransform.h', 'include/rxdock/Annotation.h',
    'include/rxdock/AnnotationHandler.h', 'include/rxdock/AromIdxSF.h',
    'include/rxdock/AtomFuncs.h', 'include/rxdock/Atom.h',
    'include/rxdock/BaseBiMolTransform.h', 'include/rxdock/BaseFileSink.h',
    'include/rxdock/BaseFileSource.h', 'include/rxdock/BaseGrid.h',
    'include/rxdock/BaseIdxSF.h', 'include/rxdock/BaseInterSF.h',
    'include/rxdock/BaseIntraSF.h', 'include/rxdock/BaseMolecularDataSource.h',
    'include/rxdock/BaseMolecularFileSink.h',
    'include/rxdock/BaseMolecularFileSource.h', 'include/rxdock/BaseObject.h',
    'include/rxdock/BaseSF.h', 'include/rxdock/BaseTransform.h',
    'include/rxdock/BaseUniMolTransform.h', 'include/rxdock/BiMolWorkSpace.h',
    'include/rxdock/Bond.h', 'include/rxdock/CavityFillSF.h',
    'include/rxdock/CavityGridSF.h', 'include/rxdock/Cavity.h',
    'include/rxdock/CellTokenIter.h', 'include/rxdock/CharmmDataSource.h',
    'include/rxdock/CharmmTypesFileSource.h',
    'include/rxdock/ChromDihedralElement.h',
    'include/rxdock/ChromDihedralRefData.h', 'include/rxdock/ChromElement.h',
    'include/rxdock/ChromFactory.h', 'include/rxdock/Chrom.h',
    'include/rxdock/ChromOccupancyElement.h',
    'include/rxdock/ChromOccupancyRefData.h',
    'include/rxdock/ChromPositionElement.h',
    'include/rxdock/ChromPositionRefData.h', 'include/rxdock/Commands.h',
    'include/rxdock/Config.h', 'include/rxdock/Constraint.h',
    'include/rxdock/ConstSF.h', 'include/rxdock/Context.h',
    'include/rxdock/Coord.h', 'include/rxdock/CrdFileSink.h',
    'include/rxdock/CrdFileSource.h', 'include/rxdock/CSVFileSink.h',
    'include/rxdock/Debug.h', 'include/rxdock/DihedralIntraSF.h',
    'include/rxdock/DihedralSF.h', 'include/rxdock/DihedralTargetSF.h',
    'include/rxdock/DirectorySource.h', 'include/rxdock/DockingError.h',
    'include/rxdock/DockingSite.h', 'include/rxdock/ElementFileSource.h',
    'include/rxdock/Error.h', 'include/rxdock/Euler.h',
    'include/rxdock/FFTGrid.h', 'include/rxdock/FileError.h',
    'include/rxdock/FilterExpression.h',
    'include/rxdock/FilterExpressionVisitor.h', 'include/rxdock/Filter.h',
    'include/rxdock/FlexAtomFactory.h', 'include/rxdock/FlexData.h',
    'include/rxdock/FlexDataVisitor.h', 'include/rxdock/GATransform.h',
    'include/rxdock/Genome.h', 'include/rxdock/InteractionGrid.h',
    'include/rxdock/InteractionTemplate.h', 'include/rxdock/LigandError.h',
    'include/rxdock/LigandFlexData.h', 'include/rxdock/LigandSiteMapper.h',
    'include/rxdock/MdlFileSink.h', 'include/rxdock/MdlFileSource.h',
    'include/rxdock/ModelError.h', 'include/rxdock/Model.h',
    'include/rxdock/ModelMutator.h', 'include/rxdock/MOEGrid.h',
    'include/rxdock/MOL2FileSource.h', 'include/rxdock/NMCriteria.h',
    'include/rxdock/NmrRestraintFileSource.h', 'include/rxdock/NmrSF.h',
    'include/rxdock/NMSimplex.h', 'include/rxdock/NMState.h',
    'include/rxdock/NoeRestraint.h', 'include/rxdock/NonBondedGrid.h',
    'include/rxdock/NonBondedHHSGrid.h', 'include/rxdock/NullTransform.h',
    'include/rxdock/Observer.h', 'include/rxdock/ParameterFileSource.h',
    'include/rxdock/ParamHandler.h', 'include/rxdock/Parser.h',
    'include/rxdock/PdbFileSource.h', 'include/rxdock/PharmaSF.h',
    'include/rxdock/Plane.h', 'include/rxdock/PMFDirSource.h',
    'include/rxdock/PMFGridSF.h', 'include/rxdock/PMF.h',
    'include/rxdock/PMFIdxSF.h', 'include/rxdock/PolarIdxSF.h',
    'include/rxdock/PolarIntraSF.h', 'include/rxdock/PolarSF.h',
    'include/rxdock/Population.h', 'include/rxdock/PrincipalAxes.h',
    'include/rxdock/PRMFactory.h', 'include/rxdock/PseudoAtom.h',
    'include/rxdock/PsfFileSink.h', 'include/rxdock/PsfFileSource.h',
    'include/rxdock/Quat.h', 'include/rxdock/Rand.h',
    'include/rxdock/RandLigTransform.h', 'include/rxdock/RandPopTransform.h',
    'include/rxdock/Rbt.h', 'include/rxdock/RealGrid.h',
    'include/rxdock/ReceptorFlexData.h', 'include/rxdock/Request.h',
    'include/rxdock/RequestHandler.h', 'include/rxdock/Resources.h',
    'include/rxdock/RotSF.h', 'include/rxdock/SAIdxSF.h',
    'include/rxdock/SATypes.h', 'include/rxdock/SetupPMFSF.h',
    'include/rxdock/SetupPolarSF.h', 'include/rxdock/SetupSASF.h',
    'include/rxdock/SFAgg.h', 'include/rxdock/SFFactory.h',
    'include/rxdock/SFRequest.h', 'include/rxdock/SimAnnTransform.h',
    'include/rxdock/SimplexCostFunction.h', 'include/rxdock/SimplexTransform.h',
    'include/rxdock/Singleton.h', 'include/rxdock/SiteMapperFactory.h',
    'include/rxdock/SiteMapper.h', 'include/rxdock/SmartPointer.h',
    'include/rxdock/SolventFlexData.h', 'include/rxdock/SphereSiteMapper.h',
    'include/rxdock/StringTokenIter.h', 'include/rxdock/Subject.h',
    'include/rxdock/TetherSF.h', 'include/rxdock/Token.h',
    'include/rxdock/TokenIter.h', 'include/rxdock/TransformAgg.h',
    'include/rxdock/TransformFactory.h', 'include/rxdock/TriposAtomType.h',
    'include/rxdock/Variant.h', 'include/rxdock/Vble.h',
    'include/rxdock/VdwGridSF.h', 'include/rxdock/VdwIdxSF.h',
    'include/rxdock/VdwIntraSF.h', 'include/rxdock/VdwSF.h',
    'include/rxdock/WorkSpace.h'),
  subdir: 'rxdock'
)
install_headers(
  files('include/rxdock/geneticprogram/Cell.h',
    'include/rxdock/geneticprogram/Command.h',
    'include/rxdock/geneticprogram/GPChromosome.h',
    'include/rxdock/geneticprogram/GPFFCHK1.h',
    'include/rxdock/geneticprogram/GPFFCrossDock.h',
    'include/rxdock/geneticprogram/GPFFGold.h',
    'include/rxdock/geneticprogram/GPFFHSP90.h',
    'include/rxdock/geneticprogram/GPFFSpike.h',
    'include/rxdock/geneticprogram/GPFitnessFunction.h',
    'include/rxdock/geneticprogram/GPGenome.h',
    'include/rxdock/geneticprogram/GPParser.h',
    'include/rxdock/geneticprogram/GPPopulation.h',
    'include/rxdock/geneticprogram/GPTypes.h'),
  subdir: 'rxdock/geneticprogram'
)

srcRbt = [
  'lib/geneticprogram/GPChromosome.cxx', 'lib/geneticprogram/GPFFCHK1.cxx',
  'lib/geneticprogram/GPFFCrossDock.cxx', 'lib/geneticprogram/GPFFGold.cxx',
  'lib/geneticprogram/GPFFHSP90.cxx', 'lib/geneticprogram/GPFFSpike.cxx',
  'lib/geneticprogram/GPFitnessFunction.cxx', 'lib/geneticprogram/GPGenome.cxx',
  'lib/geneticprogram/GPParser.cxx', 'lib/geneticprogram/GPPopulation.cxx',
  'lib/AlignTransform.cxx', 'lib/Annotation.cxx',
  'lib/AnnotationHandler.cxx', 'lib/AromIdxSF.cxx',
  'lib/Atom.cxx', 'lib/AtomFuncs.cxx',
  'lib/BaseBiMolTransform.cxx', 'lib/BaseFileSink.cxx',
  'lib/BaseFileSource.cxx', 'lib/BaseGrid.cxx',
  'lib/BaseIdxSF.cxx', 'lib/BaseInterSF.cxx',
  'lib/BaseIntraSF.cxx', 'lib/BaseMolecularFileSink.cxx',
  'lib/BaseMolecularFileSource.cxx', 'lib/BaseObject.cxx',
  'lib/BaseSF.cxx', 'lib/BaseTransform.cxx',
  'lib/BaseUniMolTransform.cxx', 'lib/BiMolWorkSpace.cxx',
  'lib/Bond.cxx', 'lib/CavityFillSF.cxx',
  'lib/CavityGridSF.cxx', 'lib/CellTokenIter.cxx',
  'lib/CharmmDataSource.cxx', 'lib/CharmmTypesFileSource.cxx',
  'lib/Chrom.cxx', 'lib/ChromDihedralElement.cxx',
  'lib/ChromDihedralRefData.cxx', 'lib/ChromElement.cxx',
  'lib/ChromFactory.cxx', 'lib/ChromOccupancyElement.cxx',
  'lib/ChromOccupancyRefData.cxx', 'lib/ChromPositionElement.cxx',
  'lib/ChromPositionRefData.cxx', 'lib/Constraint.cxx',
  'lib/ConstSF.cxx', 'lib/Context.cxx',
  'lib/CrdFileSink.cxx', 'lib/CrdFileSource.cxx',
  'lib/CSVFileSink.cxx',
  'lib/Rbt.cxx', 'lib/Debug.cxx', 'lib/DihedralIntraSF.cxx',
  'lib/DihedralSF.cxx', 'lib/DihedralTargetSF.cxx',
  'lib/DirectorySource.cxx', 'lib/DockingSite.cxx',
  'lib/ElementFileSource.cxx', 'lib/Euler.cxx',
  'lib/FFTGrid.cxx', 'lib/Filter.cxx',
  'lib/FilterExpression.cxx', 'lib/FilterExpressionVisitor.cxx',
  'lib/FlexAtomFactory.cxx', 'lib/GATransform.cxx',
  'lib/Genome.cxx', 'lib/InteractionGrid.cxx',
  'lib/LigandFlexData.cxx', 'lib/LigandSiteMapper.cxx',
  'lib/MdlFileSink.cxx', 'lib/MdlFileSource.cxx',
  'lib/Model.cxx', 'lib/ModelMutator.cxx',
  'lib/MOEGrid.cxx', 'lib/MOL2FileSource.cxx',
  'lib/NmrRestraintFileSource.cxx', 'lib/NmrSF.cxx',
  'lib/NoeRestraint.cxx', 'lib/NonBondedGrid.cxx',
  'lib/NonBondedHHSGrid.cxx', 'lib/NullTransform.cxx',
  'lib/Observer.cxx', 'lib/ParameterFileSource.cxx',
  'lib/ParamHandler.cxx', 'lib/Parser.cxx',
  'lib/PdbFileSource.cxx', 'lib/PharmaSF.cxx',
  'lib/PMF.cxx', 'lib/PMFDirSource.cxx',
  'lib/PMFGridSF.cxx', 'lib/PMFIdxSF.cxx',
  'lib/PolarIdxSF.cxx', 'lib/PolarIntraSF.cxx',
  'lib/PolarSF.cxx', 'lib/Population.cxx',
  'lib/PrincipalAxes.cxx', 'lib/PRMFactory.cxx',
  'lib/PseudoAtom.cxx', 'lib/PsfFileSink.cxx',
  'lib/PsfFileSource.cxx', 'lib/Rand.cxx',
  'lib/RandLigTransform.cxx', 'lib/RandPopTransform.cxx',
  'lib/RealGrid.cxx', 'lib/ReceptorFlexData.cxx',
  'lib/RotSF.cxx', 'lib/SAIdxSF.cxx',
  'lib/SATypes.cxx', 'lib/SetupPMFSF.cxx',
  'lib/SetupPolarSF.cxx', 'lib/SetupSASF.cxx',
  'lib/SFAgg.cxx', 'lib/SFFactory.cxx',
  'lib/SimAnnTransform.cxx', 'lib/SimplexTransform.cxx',
  'lib/SiteMapper.cxx', 'lib/SiteMapperFactory.cxx',
  'lib/SolventFlexData.cxx', 'lib/SphereSiteMapper.cxx',
  'lib/StringTokenIter.cxx', 'lib/Subject.cxx',
  'lib/TetherSF.cxx', 'lib/Token.cxx',
  'lib/TransformAgg.cxx', 'lib/TransformFactory.cxx',
  'lib/TriposAtomType.cxx', 'lib/VdwGridSF.cxx',
  'lib/VdwIdxSF.cxx', 'lib/VdwIntraSF.cxx',
  'lib/VdwSF.cxx', 'lib/WorkSpace.cxx'
]

cpp_compiler = meson.get_compiler('cpp')

eigen3_dep = dependency('eigen3', fallback : ['eigen', 'eigen_dep'])
openmp_dep = dependency('openmp', required : false)
nlohmann_json_dep = dependency('nlohmann_json', fallback : ['nlohmann_json', 'nlohmann_json_dep'])
fmt_dep = dependency('fmt', fallback : ['fmt', 'fmt_dep'])

emilk_loguru_options = ['loguru_use_fmtlib=enabled']
# FreeBSD doesn't have the required symbols to support stacktraces while the
# wrap enables them by default, so they have to be explicitely disabled
# Related pull request: https://github.com/emilk/loguru/pull/125
# Windows doesn't have libdl, which is a dependency when stactraces are enabled
if host_machine.system() == 'freebsd' or host_machine.system() == 'windows'
  emilk_loguru_options += ['loguru_stacktraces=false']
endif
emilk_loguru_dep = dependency('emilk-loguru', fallback : ['emilk-loguru', 'emilk_loguru_dep'],
                              default_options: emilk_loguru_options)

if host_machine.system() == 'windows'
  tronkko_dirent_proj = subproject('tronkko-dirent')
  tronkko_dirent_dep = tronkko_dirent_proj.get_variable('tronkko_dirent_dep')
else
  tronkko_dirent_dep = declare_dependency()
endif

# Use standard C++11 RNG on Solaris and Windows/MSVC due to PCG build failure
# Solaris issue: https://github.com/imneme/pcg-cpp/issues/42
# Windows/MSVC issue: https://github.com/imneme/pcg-cpp/issues/11
if (host_machine.system() == 'sunos' or (host_machine.system() == 'windows' and
                                         cpp_compiler.get_id() == 'msvc'))
  pcg_cpp_dep = declare_dependency()
else
  pcg_cpp_dep = dependency('pcg-cpp', fallback : ['pcg', 'pcg_cpp_dep'])
endif

library_soversion = meson.project_version().split('.')[0]
librxdock = library(
  'rxdock', srcRbt,
  dependencies : [eigen3_dep, openmp_dep, nlohmann_json_dep, pcg_cpp_dep, fmt_dep, emilk_loguru_dep, tronkko_dirent_dep],
  soversion : library_soversion,
  version : meson.project_version(),
  include_directories : incRbt, install : true
)

pkgc = import('pkgconfig')
pkgc.generate(librxdock,
    name: 'rxdock',
    version: meson.project_version(),
    description: 'A fast, versatile, and open-source library for docking ligands to proteins and nucleic acids',
)

cxxopts_dep = dependency('cxxopts', fallback : ['cxxopts', 'cxxopts_dep'])

rbdock = executable(
  'rbdock', 'tools/legacy/rbdock.cxx', link_with : librxdock,
  dependencies : [cxxopts_dep, eigen3_dep, nlohmann_json_dep, pcg_cpp_dep, fmt_dep, emilk_loguru_dep],
  include_directories : incRbt, install : true
)
rbcavity = executable(
  'rbcavity', 'tools/legacy/rbcavity.cxx', link_with : librxdock,
  dependencies : [cxxopts_dep, eigen3_dep, nlohmann_json_dep, pcg_cpp_dep],
  include_directories : incRbt, install : true
)
executable(
  'rbrms', 'tools/legacy/rbrms.cxx', link_with : librxdock,
  dependencies : [eigen3_dep, nlohmann_json_dep],
  include_directories : incRbt, install: true
)
executable(
  'rbmoegrid', 'tools/legacy/rbmoegrid.cxx', link_with : librxdock,
  dependencies : [cxxopts_dep, eigen3_dep, nlohmann_json_dep, pcg_cpp_dep],
  include_directories : incRbt, install : true
)
executable(
  'rblist', 'tools/legacy/rblist.cxx', link_with : librxdock,
  dependencies : [eigen3_dep, nlohmann_json_dep],
  include_directories : incRbt, install : true
)
executable(
  'rbcalcgrid', 'tools/legacy/rbcalcgrid.cxx', link_with : librxdock,
  dependencies : [eigen3_dep, nlohmann_json_dep, pcg_cpp_dep], include_directories : incRbt,
  install : true
)
executable(
  'rbconvgrid', 'tools/legacy/rbconvgrid.cxx', link_with : librxdock,
  dependencies : [eigen3_dep, nlohmann_json_dep], include_directories : incRbt, install : true
)

dt_smarts_dep = dependency('dt_smarts', required : false)
dt_smiles_dep = dependency('dt_smiles', required : false)
if dt_smarts_dep.found() and dt_smiles_dep.found()
  install_headers(
    files('include/rxdock/daylight/Smarts.h'),
    subdir: 'rxdock/daylight'
  )
  libRbtSmarts = library(
    'RbtSmarts', 'lib/daylight/Smarts.cxx', soversion : library_soversion,
    version : library_soversion + '.' + meson.project_name(),
    dependencies : [dt_smarts_dep, dt_smiles_dep], include_directories : incRbt
  )
  executable(
    'rbtether', 'tools/legacy/rbtether.cxx',
    dependencies : [dt_smarts_dep, dt_smiles_dep],
    link_with : [librxdock, libRbtSmarts], include_directories : incRbt
  )
  executable(
    'smart_rms', 'tools/legacy/smart_rms.cxx',
    dependencies : [dt_smarts_dep, dt_smiles_dep],
    link_with : [librxdock, libRbtSmarts], include_directories : incRbt
  )
endif

install_data([
    'data/cnx.rtf', 'data/masses.rtf', 'data/Elements.dat',
    'data/top_all22_prot_na.inp'
  ],
  install_dir : get_option('datadir') + '/' + meson.project_name().to_lower()
)

install_data('data/filters/gpfilter295.filter',
  install_dir : get_option('datadir') + '/' + meson.project_name().to_lower() +
                '/' + 'filters'
)

install_data([
    'data/pmf/CFBr.pmf', 'data/pmf/CFC0.pmf', 'data/pmf/CFC3.pmf',
    'data/pmf/CFCF.pmf', 'data/pmf/CFCL.pmf', 'data/pmf/CFCN.pmf',
    'data/pmf/CFCO.pmf', 'data/pmf/CFCP.pmf', 'data/pmf/CFCW.pmf',
    'data/pmf/CFFe.pmf', 'data/pmf/CFF_.pmf', 'data/pmf/CFHL.pmf',
    'data/pmf/CFMg.pmf', 'data/pmf/CFMn.pmf', 'data/pmf/CFN0.pmf',
    'data/pmf/CFNA.pmf', 'data/pmf/CFNC.pmf', 'data/pmf/CFND.pmf',
    'data/pmf/CFNP.pmf', 'data/pmf/CFNR.pmf', 'data/pmf/CFNS.pmf',
    'data/pmf/CFOA.pmf', 'data/pmf/CFOC.pmf', 'data/pmf/CFOD.pmf',
    'data/pmf/CFOE.pmf', 'data/pmf/CFOR.pmf', 'data/pmf/CFOS.pmf',
    'data/pmf/CFP_.pmf', 'data/pmf/CFSA.pmf', 'data/pmf/CFSD.pmf',
    'data/pmf/CFV_.pmf', 'data/pmf/CFZn.pmf', 'data/pmf/CNBr.pmf',
    'data/pmf/CNC0.pmf', 'data/pmf/CNC3.pmf', 'data/pmf/CNCF.pmf',
    'data/pmf/CNCL.pmf', 'data/pmf/CNCN.pmf', 'data/pmf/CNCO.pmf',
    'data/pmf/CNCP.pmf', 'data/pmf/CNCW.pmf', 'data/pmf/CNFe.pmf',
    'data/pmf/CNF_.pmf', 'data/pmf/CNHL.pmf', 'data/pmf/CNMg.pmf',
    'data/pmf/CNMn.pmf', 'data/pmf/CNN0.pmf', 'data/pmf/CNNA.pmf',
    'data/pmf/CNNC.pmf', 'data/pmf/CNND.pmf', 'data/pmf/CNNP.pmf',
    'data/pmf/CNNR.pmf', 'data/pmf/CNNS.pmf', 'data/pmf/CNOA.pmf',
    'data/pmf/CNOC.pmf', 'data/pmf/CNOD.pmf', 'data/pmf/CNOE.pmf',
    'data/pmf/CNOR.pmf', 'data/pmf/CNOS.pmf', 'data/pmf/CNP_.pmf',
    'data/pmf/CNSA.pmf', 'data/pmf/CNSD.pmf', 'data/pmf/CNV_.pmf',
    'data/pmf/CNZn.pmf', 'data/pmf/COBr.pmf', 'data/pmf/COC0.pmf',
    'data/pmf/COC3.pmf', 'data/pmf/COCF.pmf', 'data/pmf/COCL.pmf',
    'data/pmf/COCN.pmf', 'data/pmf/COCO.pmf', 'data/pmf/COCP.pmf',
    'data/pmf/COCW.pmf', 'data/pmf/COFe.pmf', 'data/pmf/COF_.pmf',
    'data/pmf/COHL.pmf', 'data/pmf/COMg.pmf', 'data/pmf/COMn.pmf',
    'data/pmf/CON0.pmf', 'data/pmf/CONA.pmf', 'data/pmf/CONC.pmf',
    'data/pmf/COND.pmf', 'data/pmf/CONP.pmf', 'data/pmf/CONR.pmf',
    'data/pmf/CONS.pmf', 'data/pmf/COOA.pmf', 'data/pmf/COOC.pmf',
    'data/pmf/COOD.pmf', 'data/pmf/COOE.pmf', 'data/pmf/COOR.pmf',
    'data/pmf/COOS.pmf', 'data/pmf/COP_.pmf', 'data/pmf/COSA.pmf',
    'data/pmf/COSD.pmf', 'data/pmf/COV_.pmf', 'data/pmf/COZn.pmf',
    'data/pmf/CPBr.pmf', 'data/pmf/CPC0.pmf', 'data/pmf/CPC3.pmf',
    'data/pmf/CPCF.pmf', 'data/pmf/CPCL.pmf', 'data/pmf/CPCN.pmf',
    'data/pmf/CPCO.pmf', 'data/pmf/CPCP.pmf', 'data/pmf/CPCW.pmf',
    'data/pmf/CPFe.pmf', 'data/pmf/CPF_.pmf', 'data/pmf/CPHL.pmf',
    'data/pmf/CPMg.pmf', 'data/pmf/CPMn.pmf', 'data/pmf/CPN0.pmf',
    'data/pmf/CPNA.pmf', 'data/pmf/CPNC.pmf', 'data/pmf/CPND.pmf',
    'data/pmf/CPNP.pmf', 'data/pmf/CPNR.pmf', 'data/pmf/CPNS.pmf',
    'data/pmf/CPOA.pmf', 'data/pmf/CPOC.pmf', 'data/pmf/CPOD.pmf',
    'data/pmf/CPOE.pmf', 'data/pmf/CPOR.pmf', 'data/pmf/CPOS.pmf',
    'data/pmf/CPP_.pmf', 'data/pmf/CPSA.pmf', 'data/pmf/CPSD.pmf',
    'data/pmf/CPV_.pmf', 'data/pmf/CPZn.pmf', 'data/pmf/HHBr.pmf',
    'data/pmf/HHC0.pmf', 'data/pmf/HHC3.pmf', 'data/pmf/HHCF.pmf',
    'data/pmf/HHCL.pmf', 'data/pmf/HHCN.pmf', 'data/pmf/HHCO.pmf',
    'data/pmf/HHCP.pmf', 'data/pmf/HHCW.pmf', 'data/pmf/HHFe.pmf',
    'data/pmf/HHF_.pmf', 'data/pmf/HHHL.pmf', 'data/pmf/HHMg.pmf',
    'data/pmf/HHMn.pmf', 'data/pmf/HHN0.pmf', 'data/pmf/HHNA.pmf',
    'data/pmf/HHNC.pmf', 'data/pmf/HHND.pmf', 'data/pmf/HHNP.pmf',
    'data/pmf/HHNR.pmf', 'data/pmf/HHNS.pmf', 'data/pmf/HHOA.pmf',
    'data/pmf/HHOC.pmf', 'data/pmf/HHOD.pmf', 'data/pmf/HHOE.pmf',
    'data/pmf/HHOR.pmf', 'data/pmf/HHOS.pmf', 'data/pmf/HHP_.pmf',
    'data/pmf/HHSA.pmf', 'data/pmf/HHSD.pmf', 'data/pmf/HHV_.pmf',
    'data/pmf/HHZn.pmf', 'data/pmf/NCBr.pmf', 'data/pmf/NCC0.pmf',
    'data/pmf/NCC3.pmf', 'data/pmf/NCCF.pmf', 'data/pmf/NCCL.pmf',
    'data/pmf/NCCN.pmf', 'data/pmf/NCCO.pmf', 'data/pmf/NCCP.pmf',
    'data/pmf/NCCW.pmf', 'data/pmf/NCFe.pmf', 'data/pmf/NCF_.pmf',
    'data/pmf/NCHL.pmf', 'data/pmf/NCMg.pmf', 'data/pmf/NCMn.pmf',
    'data/pmf/NCN0.pmf', 'data/pmf/NCNA.pmf', 'data/pmf/NCNC.pmf',
    'data/pmf/NCND.pmf', 'data/pmf/NCNP.pmf', 'data/pmf/NCNR.pmf',
    'data/pmf/NCNS.pmf', 'data/pmf/NCOA.pmf', 'data/pmf/NCOC.pmf',
    'data/pmf/NCOD.pmf', 'data/pmf/NCOE.pmf', 'data/pmf/NCOR.pmf',
    'data/pmf/NCOS.pmf', 'data/pmf/NCP_.pmf', 'data/pmf/NCSA.pmf',
    'data/pmf/NCSD.pmf', 'data/pmf/NCV_.pmf', 'data/pmf/NCZn.pmf',
    'data/pmf/NDBr.pmf', 'data/pmf/NDC0.pmf', 'data/pmf/NDC3.pmf',
    'data/pmf/NDCF.pmf', 'data/pmf/NDCL.pmf', 'data/pmf/NDCN.pmf',
    'data/pmf/NDCO.pmf', 'data/pmf/NDCP.pmf', 'data/pmf/NDCW.pmf',
    'data/pmf/NDFe.pmf', 'data/pmf/NDF_.pmf', 'data/pmf/NDHL.pmf',
    'data/pmf/NDMg.pmf', 'data/pmf/NDMn.pmf', 'data/pmf/NDN0.pmf',
    'data/pmf/NDNA.pmf', 'data/pmf/NDNC.pmf', 'data/pmf/NDND.pmf',
    'data/pmf/NDNP.pmf', 'data/pmf/NDNR.pmf', 'data/pmf/NDNS.pmf',
    'data/pmf/NDOA.pmf', 'data/pmf/NDOC.pmf', 'data/pmf/NDOD.pmf',
    'data/pmf/NDOE.pmf', 'data/pmf/NDOR.pmf', 'data/pmf/NDOS.pmf',
    'data/pmf/NDP_.pmf', 'data/pmf/NDSA.pmf', 'data/pmf/NDSD.pmf',
    'data/pmf/NDV_.pmf', 'data/pmf/NDZn.pmf', 'data/pmf/NRBr.pmf',
    'data/pmf/NRC0.pmf', 'data/pmf/NRC3.pmf', 'data/pmf/NRCF.pmf',
    'data/pmf/NRCL.pmf', 'data/pmf/NRCN.pmf', 'data/pmf/NRCO.pmf',
    'data/pmf/NRCP.pmf', 'data/pmf/NRCW.pmf', 'data/pmf/NRFe.pmf',
    'data/pmf/NRF_.pmf', 'data/pmf/NRHL.pmf', 'data/pmf/NRMg.pmf',
    'data/pmf/NRMn.pmf', 'data/pmf/NRN0.pmf', 'data/pmf/NRNA.pmf',
    'data/pmf/NRNC.pmf', 'data/pmf/NRND.pmf', 'data/pmf/NRNP.pmf',
    'data/pmf/NRNR.pmf', 'data/pmf/NRNS.pmf', 'data/pmf/NROA.pmf',
    'data/pmf/NROC.pmf', 'data/pmf/NROD.pmf', 'data/pmf/NROE.pmf',
    'data/pmf/NROR.pmf', 'data/pmf/NROS.pmf', 'data/pmf/NRP_.pmf',
    'data/pmf/NRSA.pmf', 'data/pmf/NRSD.pmf', 'data/pmf/NRV_.pmf',
    'data/pmf/NRZn.pmf', 'data/pmf/OABr.pmf', 'data/pmf/OAC0.pmf',
    'data/pmf/OAC3.pmf', 'data/pmf/OACF.pmf', 'data/pmf/OACL.pmf',
    'data/pmf/OACN.pmf', 'data/pmf/OACO.pmf', 'data/pmf/OACP.pmf',
    'data/pmf/OACW.pmf', 'data/pmf/OAFe.pmf', 'data/pmf/OAF_.pmf',
    'data/pmf/OAHL.pmf', 'data/pmf/OAMg.pmf', 'data/pmf/OAMn.pmf',
    'data/pmf/OAN0.pmf', 'data/pmf/OANA.pmf', 'data/pmf/OANC.pmf',
    'data/pmf/OAND.pmf', 'data/pmf/OANP.pmf', 'data/pmf/OANR.pmf',
    'data/pmf/OANS.pmf', 'data/pmf/OAOA.pmf', 'data/pmf/OAOC.pmf',
    'data/pmf/OAOD.pmf', 'data/pmf/OAOE.pmf', 'data/pmf/OAOR.pmf',
    'data/pmf/OAOS.pmf', 'data/pmf/OAP_.pmf', 'data/pmf/OASA.pmf',
    'data/pmf/OASD.pmf', 'data/pmf/OAV_.pmf', 'data/pmf/OAZn.pmf',
    'data/pmf/OCBr.pmf', 'data/pmf/OCC0.pmf', 'data/pmf/OCC3.pmf',
    'data/pmf/OCCF.pmf', 'data/pmf/OCCL.pmf', 'data/pmf/OCCN.pmf',
    'data/pmf/OCCO.pmf', 'data/pmf/OCCP.pmf', 'data/pmf/OCCW.pmf',
    'data/pmf/OCFe.pmf', 'data/pmf/OCF_.pmf', 'data/pmf/OCHL.pmf',
    'data/pmf/OCMg.pmf', 'data/pmf/OCMn.pmf', 'data/pmf/OCN0.pmf',
    'data/pmf/OCNA.pmf', 'data/pmf/OCNC.pmf', 'data/pmf/OCND.pmf',
    'data/pmf/OCNP.pmf', 'data/pmf/OCNR.pmf', 'data/pmf/OCNS.pmf',
    'data/pmf/OCOA.pmf', 'data/pmf/OCOC.pmf', 'data/pmf/OCOD.pmf',
    'data/pmf/OCOE.pmf', 'data/pmf/OCOR.pmf', 'data/pmf/OCOS.pmf',
    'data/pmf/OCP_.pmf', 'data/pmf/OCSA.pmf', 'data/pmf/OCSD.pmf',
    'data/pmf/OCV_.pmf', 'data/pmf/OCZn.pmf', 'data/pmf/ODBr.pmf',
    'data/pmf/ODC0.pmf', 'data/pmf/ODC3.pmf', 'data/pmf/ODCF.pmf',
    'data/pmf/ODCL.pmf', 'data/pmf/ODCN.pmf', 'data/pmf/ODCO.pmf',
    'data/pmf/ODCP.pmf', 'data/pmf/ODCW.pmf', 'data/pmf/ODFe.pmf',
    'data/pmf/ODF_.pmf', 'data/pmf/ODHL.pmf', 'data/pmf/ODMg.pmf',
    'data/pmf/ODMn.pmf', 'data/pmf/ODN0.pmf', 'data/pmf/ODNA.pmf',
    'data/pmf/ODNC.pmf', 'data/pmf/ODND.pmf', 'data/pmf/ODNP.pmf',
    'data/pmf/ODNR.pmf', 'data/pmf/ODNS.pmf', 'data/pmf/ODOA.pmf',
    'data/pmf/ODOC.pmf', 'data/pmf/ODOD.pmf', 'data/pmf/ODOE.pmf',
    'data/pmf/ODOR.pmf', 'data/pmf/ODOS.pmf', 'data/pmf/ODP_.pmf',
    'data/pmf/ODSA.pmf', 'data/pmf/ODSD.pmf', 'data/pmf/ODV_.pmf',
    'data/pmf/ODZn.pmf', 'data/pmf/OWBr.pmf', 'data/pmf/OWC0.pmf',
    'data/pmf/OWC3.pmf', 'data/pmf/OWCF.pmf', 'data/pmf/OWCL.pmf',
    'data/pmf/OWCN.pmf', 'data/pmf/OWCO.pmf', 'data/pmf/OWCP.pmf',
    'data/pmf/OWCW.pmf', 'data/pmf/OWFe.pmf', 'data/pmf/OWF_.pmf',
    'data/pmf/OWHL.pmf', 'data/pmf/OWMg.pmf', 'data/pmf/OWMn.pmf',
    'data/pmf/OWN0.pmf', 'data/pmf/OWNA.pmf', 'data/pmf/OWNC.pmf',
    'data/pmf/OWND.pmf', 'data/pmf/OWNP.pmf', 'data/pmf/OWNR.pmf',
    'data/pmf/OWNS.pmf', 'data/pmf/OWOA.pmf', 'data/pmf/OWOC.pmf',
    'data/pmf/OWOD.pmf', 'data/pmf/OWOE.pmf', 'data/pmf/OWOR.pmf',
    'data/pmf/OWOS.pmf', 'data/pmf/OWP_.pmf', 'data/pmf/OWSA.pmf',
    'data/pmf/OWSD.pmf', 'data/pmf/OWV_.pmf', 'data/pmf/OWZn.pmf',
    'data/pmf/SABr.pmf', 'data/pmf/SAC0.pmf', 'data/pmf/SAC3.pmf',
    'data/pmf/SACF.pmf', 'data/pmf/SACL.pmf', 'data/pmf/SACN.pmf',
    'data/pmf/SACO.pmf', 'data/pmf/SACP.pmf', 'data/pmf/SACW.pmf',
    'data/pmf/SAFe.pmf', 'data/pmf/SAF_.pmf', 'data/pmf/SAHL.pmf',
    'data/pmf/SAMg.pmf', 'data/pmf/SAMn.pmf', 'data/pmf/SAN0.pmf',
    'data/pmf/SANA.pmf', 'data/pmf/SANC.pmf', 'data/pmf/SAND.pmf',
    'data/pmf/SANP.pmf', 'data/pmf/SANR.pmf', 'data/pmf/SANS.pmf',
    'data/pmf/SAOA.pmf', 'data/pmf/SAOC.pmf', 'data/pmf/SAOD.pmf',
    'data/pmf/SAOE.pmf', 'data/pmf/SAOR.pmf', 'data/pmf/SAOS.pmf',
    'data/pmf/SAP_.pmf', 'data/pmf/SASA.pmf', 'data/pmf/SASD.pmf',
    'data/pmf/SAV_.pmf', 'data/pmf/SAZn.pmf', 'data/pmf/SDBr.pmf',
    'data/pmf/SDC0.pmf', 'data/pmf/SDC3.pmf', 'data/pmf/SDCF.pmf',
    'data/pmf/SDCL.pmf', 'data/pmf/SDCN.pmf', 'data/pmf/SDCO.pmf',
    'data/pmf/SDCP.pmf', 'data/pmf/SDCW.pmf', 'data/pmf/SDFe.pmf',
    'data/pmf/SDF_.pmf', 'data/pmf/SDHL.pmf', 'data/pmf/SDMg.pmf',
    'data/pmf/SDMn.pmf', 'data/pmf/SDN0.pmf', 'data/pmf/SDNA.pmf',
    'data/pmf/SDNC.pmf', 'data/pmf/SDND.pmf', 'data/pmf/SDNP.pmf',
    'data/pmf/SDNR.pmf', 'data/pmf/SDNS.pmf', 'data/pmf/SDOA.pmf',
    'data/pmf/SDOC.pmf', 'data/pmf/SDOD.pmf', 'data/pmf/SDOE.pmf',
    'data/pmf/SDOR.pmf', 'data/pmf/SDOS.pmf', 'data/pmf/SDP_.pmf',
    'data/pmf/SDSA.pmf', 'data/pmf/SDSD.pmf', 'data/pmf/SDV_.pmf',
    'data/pmf/SDZn.pmf'
  ],
  install_dir : get_option('datadir') + '/' + meson.project_name().to_lower() +
                '/' + 'pmf'
)

install_data([
    'data/scripts/dock_grid.prm', 'data/scripts/dock.prm',
    'data/scripts/dock_solv_grid.prm', 'data/scripts/dock_solv.prm',
    'data/scripts/minimise.prm', 'data/scripts/minimise_solv.prm',
    'data/scripts/score_pmf.prm', 'data/scripts/score.prm',
    'data/scripts/score_solv.prm'
  ],
  install_dir : get_option('datadir') + '/' + meson.project_name().to_lower() +
                '/' + 'scripts'
)

install_data([
    'data/sf/calcgrid_vdw1.prm', 'data/sf/calcgrid_vdw5.prm',
    'data/sf/InterGridSF.prm', 'data/sf/InterIdxSF.prm',
    'data/sf/IntraSF.prm', 'data/sf/IonicAtoms.prm',
    'data/sf/PMFIdxSF.prm', 'data/sf/SolvGridSF.prm',
    'data/sf/SolvIdxSF.prm', 'data/sf/TargetSF.prm',
    'data/sf/solvation_asp.prm', 'data/sf/Tripos52_dihedrals.prm',
    'data/sf/Tripos52_vdw.prm'
  ],
  install_dir : get_option('datadir') + '/' + meson.project_name().to_lower() +
                '/' + 'sf'
)

tests_opt = get_option('tests')
if tests_opt
  gtest_dep = dependency(
    'gtest', fallback : ['gtest', 'gtest_dep'], required : false
  )
  if gtest_dep.found()
    incTest = include_directories('tests')
    srcTest = [
      'tests/Main.cxx', 'tests/OccupancyTest.cxx',
      'tests/ChromTest.cxx', 'tests/SearchTest.cxx'
    ]
    unit_test = executable(
      'unit-test', srcTest,
      dependencies : [pcg_cpp_dep, gtest_dep, eigen3_dep, nlohmann_json_dep],
      link_with : librxdock, include_directories : [incTest, incRbt]
    )
    rbtHomePath = '/tests/data'
    test(
      'rbcavity-1koc', rbcavity, args : [
        '-r', meson.current_source_dir() + rbtHomePath + '/1koc.prm',
        '-W'
      ],
      env : [
        'RBT_ROOT=' + meson.current_source_dir(),
        'RBT_HOME=' + meson.current_source_dir() + rbtHomePath
      ],
      timeout : 900
    )
    test(
      'rbcavity-1yet', rbcavity, args : [
        '-r', meson.current_source_dir() + rbtHomePath + '/1YET.prm',
        '-W'
      ],
      env : [
        'RBT_ROOT=' + meson.current_source_dir(),
        'RBT_HOME=' + meson.current_source_dir() + rbtHomePath
      ],
      timeout : 900
    )
    test(
      'rbcavity-1yet-test', rbcavity, args : [
        '-r', meson.current_source_dir() + rbtHomePath + '/1YET_test.prm',
        '-W'
      ],
      env : [
        'RBT_ROOT=' + meson.current_source_dir(),
        'RBT_HOME=' + meson.current_source_dir() + rbtHomePath
      ],
      timeout : 900
    )
    test(
      'unit-test', unit_test,
      env : [
        'RBT_ROOT=' + meson.current_source_dir(),
        'RBT_HOME=' + meson.current_source_dir() + rbtHomePath
      ],
      is_parallel : false,
      timeout : 900
    )
    test(
      'rbdock-1yet-test', rbdock, args : [
        '-r', meson.current_source_dir() + rbtHomePath + '/1YET_test.prm',
        '-i', meson.current_source_dir() + rbtHomePath + '/1YET_c.sd',
        '-p', 'dock.prm', '-n', '1', '-s', '48151623',
        '-o', meson.current_source_dir() + rbtHomePath + '/1YET_test_out'
      ],
      env : [
        'RBT_ROOT=' + meson.current_source_dir(),
        'RBT_HOME=' + meson.current_source_dir() + rbtHomePath
      ],
      is_parallel : false,
      timeout : 900
    )
  endif
endif
