project('rDock', 'cpp', version : '2013.1', default_options : ['cpp_std=c++11'])

incRbt = include_directories('include', 'include/GP')
incRbtAll = include_directories(
  'include', 'include/GP', 'import/simplex/include', 'import/tnt/include'
)
srcRbt = [
  'src/GP/RbtGPChromosome.cxx', 'src/GP/RbtGPFFCHK1.cxx',
  'src/GP/RbtGPFFCrossDock.cxx', 'src/GP/RbtGPFFGold.cxx',
  'src/GP/RbtGPFFHSP90.cxx', 'src/GP/RbtGPFFSpike.cxx',
  'src/GP/RbtGPFitnessFunction.cxx', 'src/GP/RbtGPGenome.cxx',
  'src/GP/RbtGPParser.cxx', 'src/GP/RbtGPPopulation.cxx',
  'src/lib/RbtAlignTransform.cxx', 'src/lib/RbtAnnotation.cxx',
  'src/lib/RbtAnnotationHandler.cxx', 'src/lib/RbtAromIdxSF.cxx',
  'src/lib/RbtAtom.cxx', 'src/lib/RbtAtomFuncs.cxx',
  'src/lib/RbtBaseBiMolTransform.cxx', 'src/lib/RbtBaseFileSink.cxx',
  'src/lib/RbtBaseFileSource.cxx', 'src/lib/RbtBaseGrid.cxx',
  'src/lib/RbtBaseIdxSF.cxx', 'src/lib/RbtBaseInterSF.cxx',
  'src/lib/RbtBaseIntraSF.cxx', 'src/lib/RbtBaseMolecularFileSink.cxx',
  'src/lib/RbtBaseMolecularFileSource.cxx', 'src/lib/RbtBaseObject.cxx',
  'src/lib/RbtBaseSF.cxx', 'src/lib/RbtBaseTransform.cxx',
  'src/lib/RbtBaseUniMolTransform.cxx', 'src/lib/RbtBiMolWorkSpace.cxx',
  'src/lib/RbtBond.cxx', 'src/lib/RbtCavityFillSF.cxx',
  'src/lib/RbtCavityGridSF.cxx', 'src/lib/RbtCellTokenIter.cxx',
  'src/lib/RbtCharmmDataSource.cxx', 'src/lib/RbtCharmmTypesFileSource.cxx',
  'src/lib/RbtChrom.cxx', 'src/lib/RbtChromDihedralElement.cxx',
  'src/lib/RbtChromDihedralRefData.cxx', 'src/lib/RbtChromElement.cxx',
  'src/lib/RbtChromFactory.cxx', 'src/lib/RbtChromOccupancyElement.cxx',
  'src/lib/RbtChromOccupancyRefData.cxx', 'src/lib/RbtChromPositionElement.cxx',
  'src/lib/RbtChromPositionRefData.cxx', 'src/lib/RbtConstraint.cxx',
  'src/lib/RbtConstSF.cxx', 'src/lib/RbtContext.cxx',
  'src/lib/RbtCrdFileSink.cxx', 'src/lib/RbtCrdFileSource.cxx',
  'src/lib/Rbt.cxx', 'src/lib/RbtDebug.cxx', 'src/lib/RbtDihedralIntraSF.cxx',
  'src/lib/RbtDihedralSF.cxx', 'src/lib/RbtDihedralTargetSF.cxx',
  'src/lib/RbtDirectorySource.cxx', 'src/lib/RbtDockingSite.cxx',
  'src/lib/RbtElementFileSource.cxx', 'src/lib/RbtEuler.cxx',
  'src/lib/RbtFFTGrid.cxx', 'src/lib/RbtFilter.cxx',
  'src/lib/RbtFilterExpression.cxx', 'src/lib/RbtFilterExpressionVisitor.cxx',
  'src/lib/RbtFlexAtomFactory.cxx', 'src/lib/RbtGATransform.cxx',
  'src/lib/RbtGenome.cxx', 'src/lib/RbtInteractionGrid.cxx',
  'src/lib/RbtLigandFlexData.cxx', 'src/lib/RbtLigandSiteMapper.cxx',
  'src/lib/RbtMdlFileSink.cxx', 'src/lib/RbtMdlFileSource.cxx',
  'src/lib/RbtModel.cxx', 'src/lib/RbtModelMutator.cxx',
  'src/lib/RbtMOEGrid.cxx', 'src/lib/RbtMOL2FileSource.cxx',
  'src/lib/RbtNmrRestraintFileSource.cxx', 'src/lib/RbtNmrSF.cxx',
  'src/lib/RbtNoeRestraint.cxx', 'src/lib/RbtNonBondedGrid.cxx',
  'src/lib/RbtNonBondedHHSGrid.cxx', 'src/lib/RbtNullTransform.cxx',
  'src/lib/RbtObserver.cxx', 'src/lib/RbtParameterFileSource.cxx',
  'src/lib/RbtParamHandler.cxx', 'src/lib/RbtParser.cxx',
  'src/lib/RbtPdbFileSource.cxx', 'src/lib/RbtPharmaSF.cxx',
  'src/lib/RbtPMF.cxx', 'src/lib/RbtPMFDirSource.cxx',
  'src/lib/RbtPMFGridSF.cxx', 'src/lib/RbtPMFIdxSF.cxx',
  'src/lib/RbtPolarIdxSF.cxx', 'src/lib/RbtPolarIntraSF.cxx',
  'src/lib/RbtPolarSF.cxx', 'src/lib/RbtPopulation.cxx',
  'src/lib/RbtPrincipalAxes.cxx', 'src/lib/RbtPRMFactory.cxx',
  'src/lib/RbtPseudoAtom.cxx', 'src/lib/RbtPsfFileSink.cxx',
  'src/lib/RbtPsfFileSource.cxx', 'src/lib/RbtRand.cxx',
  'src/lib/RbtRandLigTransform.cxx', 'src/lib/RbtRandPopTransform.cxx',
  'src/lib/RbtRealGrid.cxx', 'src/lib/RbtReceptorFlexData.cxx',
  'src/lib/RbtRotSF.cxx', 'src/lib/RbtSAIdxSF.cxx',
  'src/lib/RbtSATypes.cxx', 'src/lib/RbtSetupPMFSF.cxx',
  'src/lib/RbtSetupPolarSF.cxx', 'src/lib/RbtSetupSASF.cxx',
  'src/lib/RbtSFAgg.cxx', 'src/lib/RbtSFFactory.cxx',
  'src/lib/RbtSimAnnTransform.cxx', 'src/lib/RbtSimplexTransform.cxx',
  'src/lib/RbtSiteMapper.cxx', 'src/lib/RbtSiteMapperFactory.cxx',
  'src/lib/RbtSolventFlexData.cxx', 'src/lib/RbtSphereSiteMapper.cxx',
  'src/lib/RbtStringTokenIter.cxx', 'src/lib/RbtSubject.cxx',
  'src/lib/RbtTetherSF.cxx', 'src/lib/RbtToken.cxx',
  'src/lib/RbtTransformAgg.cxx', 'src/lib/RbtTransformFactory.cxx',
  'src/lib/RbtTriposAtomType.cxx', 'src/lib/RbtVdwGridSF.cxx',
  'src/lib/RbtVdwIdxSF.cxx', 'src/lib/RbtVdwIntraSF.cxx',
  'src/lib/RbtVdwSF.cxx', 'src/lib/RbtWorkSpace.cxx',
  'import/simplex/src/NMSearch.cxx'
]

pcg_cpp_dep = dependency('pcg-cpp', fallback : ['pcg', 'pcg_cpp_dep'])

library_soversion = '0'
libRbt = library(
  'Rbt', srcRbt, dependencies : pcg_cpp_dep, soversion : library_soversion,
  version : library_soversion + '.' + meson.project_version(),
  include_directories : incRbtAll
)

rbdock = executable(
  'rbdock', 'src/exe/rbdock.cxx', link_with : libRbt,
  dependencies : pcg_cpp_dep, include_directories : incRbt
)
rbcavity = executable(
  'rbcavity', 'src/exe/rbcavity.cxx', link_with : libRbt,
  dependencies : pcg_cpp_dep, include_directories : incRbt
)
executable(
  'rbrms', 'src/exe/rbrms.cxx', link_with : libRbt, include_directories : incRbt
)
executable(
  'rbmoegrid', 'src/exe/rbmoegrid.cxx', link_with : libRbt,
  dependencies : pcg_cpp_dep, include_directories : incRbt
)
executable(
  'rblist', 'src/exe/rblist.cxx', link_with : libRbt,
  include_directories : incRbt
)
executable(
  'rbcalcgrid', 'src/exe/rbcalcgrid.cxx', link_with : libRbt,
  dependencies : pcg_cpp_dep, include_directories : incRbt
)
executable(
  'rbconvgrid', 'src/exe/rbconvgrid.cxx', link_with : libRbt,
  include_directories : incRbt
)

dt_smarts_dep = dependency('dt_smarts', required : false)
dt_smiles_dep = dependency('dt_smiles', required : false)
if dt_smarts_dep.found() and dt_smiles_dep.found()
  libRbtSmarts = library(
    'RbtSmarts', 'src/daylight/RbtSmarts.cxx', soversion : library_soversion,
    version : library_soversion + '.' + meson.project_version(),
    dependencies : [dt_smarts_dep, dt_smiles_dep], include_directories : incRbt
  )
  executable(
    'rbtether', 'src/exe/rbtether.cxx',
    dependencies : [dt_smarts_dep, dt_smiles_dep],
    link_with : [libRbt, libRbtSmarts], include_directories : incRbt
  )
  executable(
    'smart_rms', 'src/exe/smart_rms.cxx',
    dependencies : [dt_smarts_dep, dt_smiles_dep],
    link_with : [libRbt, libRbtSmarts], include_directories : incRbt
  )
endif

tests_opt = get_option('tests')
if tests_opt
  gtest_dep = dependency(
    'gtest', fallback : ['gtest', 'gtest_dep'], required : false
  )
  if gtest_dep.found()
    incTest = include_directories('build/test')
    srcTest = [
      'build/test/Main.cxx', 'build/test/OccupancyTest.cxx',
      'build/test/RbtChromTest.cxx', 'build/test/SearchTest.cxx'
    ]
    unit_test = executable(
      'unit-test', srcTest, dependencies : [pcg_cpp_dep, gtest_dep],
      link_with : libRbt, include_directories : [incTest, incRbt]
    )
    rbtHomePath = '/build/test/RBT_HOME'
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
