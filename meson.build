project('RxDock', 'cpp', version : '0.1.0', default_options : ['cpp_std=c++11'])

incRbt = include_directories('include', 'include/GP')
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
  'src/lib/RbtVdwSF.cxx', 'src/lib/RbtWorkSpace.cxx'
]

cpp_compiler = meson.get_compiler('cpp')

eigen3_dep = dependency('eigen3', fallback : ['eigen', 'eigen_dep'])

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
  dependencies : [eigen3_dep, pcg_cpp_dep, tronkko_dirent_dep],
  soversion : library_soversion,
  version : meson.project_version(),
  include_directories : incRbt, install : true
)

cxxopts_dep = dependency('cxxopts', fallback : ['cxxopts', 'cxxopts_dep'])

rbdock = executable(
  'rbdock', 'src/exe/rbdock.cxx', link_with : librxdock,
  dependencies : [cxxopts_dep, eigen3_dep, pcg_cpp_dep],
  include_directories : incRbt, install : true
)
rbcavity = executable(
  'rbcavity', 'src/exe/rbcavity.cxx', link_with : librxdock,
  dependencies : [cxxopts_dep, eigen3_dep, pcg_cpp_dep],
  include_directories : incRbt, install : true
)
executable(
  'rbrms', 'src/exe/rbrms.cxx', link_with : librxdock,
  include_directories : incRbt, install: true
)
executable(
  'rbmoegrid', 'src/exe/rbmoegrid.cxx', link_with : librxdock,
  dependencies : [cxxopts_dep, eigen3_dep, pcg_cpp_dep],
  include_directories : incRbt, install : true
)
executable(
  'rblist', 'src/exe/rblist.cxx', link_with : librxdock,
  include_directories : incRbt, install : true
)
executable(
  'rbcalcgrid', 'src/exe/rbcalcgrid.cxx', link_with : librxdock,
  dependencies : [eigen3_dep, pcg_cpp_dep], include_directories : incRbt,
  install : true
)
executable(
  'rbconvgrid', 'src/exe/rbconvgrid.cxx', link_with : librxdock,
  dependencies : eigen3_dep, include_directories : incRbt, install : true
)

dt_smarts_dep = dependency('dt_smarts', required : false)
dt_smiles_dep = dependency('dt_smiles', required : false)
if dt_smarts_dep.found() and dt_smiles_dep.found()
  libRbtSmarts = library(
    'RbtSmarts', 'src/daylight/RbtSmarts.cxx', soversion : library_soversion,
    version : library_soversion + '.' + meson.project_name(),
    dependencies : [dt_smarts_dep, dt_smiles_dep], include_directories : incRbt
  )
  executable(
    'rbtether', 'src/exe/rbtether.cxx',
    dependencies : [dt_smarts_dep, dt_smiles_dep],
    link_with : [librxdock, libRbtSmarts], include_directories : incRbt
  )
  executable(
    'smart_rms', 'src/exe/smart_rms.cxx',
    dependencies : [dt_smarts_dep, dt_smiles_dep],
    link_with : [librxdock, libRbtSmarts], include_directories : incRbt
  )
endif

install_data([
    'data/cnx.rtf', 'data/masses.rtf', 'data/RbtElements.dat',
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
    'data/sf/RbtInterGridSF.prm', 'data/sf/RbtInterIdxSF.prm',
    'data/sf/RbtIntraSF.prm', 'data/sf/RbtIonicAtoms.prm',
    'data/sf/RbtPMFIdxSF.prm', 'data/sf/RbtSolvGridSF.prm',
    'data/sf/RbtSolvIdxSF.prm', 'data/sf/RbtTargetSF.prm',
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
    incTest = include_directories('build/test')
    srcTest = [
      'build/test/Main.cxx', 'build/test/OccupancyTest.cxx',
      'build/test/RbtChromTest.cxx', 'build/test/SearchTest.cxx'
    ]
    unit_test = executable(
      'unit-test', srcTest,
      dependencies : [pcg_cpp_dep, gtest_dep, eigen3_dep],
      link_with : librxdock, include_directories : [incTest, incRbt]
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
