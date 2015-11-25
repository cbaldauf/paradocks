#include "ParadocksConfig.hpp"

#include <fstream>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>

XERCES_CPP_NAMESPACE_USE

namespace paradocks
{

string XMLString2string(const XMLCh* toTranscode)
{
	char* tmp=XMLString::transcode(toTranscode);
	string result(tmp);
	XMLString::release(&tmp);
	return result;
}

ParadocksConfig::ParadocksConfig(string configfile)
{
	fstream f(configfile.c_str());
	if (!f.good())
	{
		string emsg("Configuration file \""+configfile+"\" could not be opened for reading!");
		throw ParadocksConfigError(emsg);
	}
	// initialize XML
	try
	{
		XMLPlatformUtils::Initialize();
	}
	catch (const XMLException& e)
	{
		string emsg("Error during XML initialization! :\n"+XMLString2string(e.getMessage()));
		throw ParadocksConfigError(emsg);
	}
	
	// parse config file
	XercesDOMParser* ConfigFileParser=new XercesDOMParser;
	DOMDocument* xmlDoc;
	try
	{
		ConfigFileParser->parse(configfile.c_str());
		// Configure DOM parser.
		ConfigFileParser->setValidationScheme(XercesDOMParser::Val_Never);
		ConfigFileParser->setDoNamespaces(false);
		ConfigFileParser->setDoSchema(false);
		ConfigFileParser->setLoadExternalDTD(false);
		xmlDoc = ConfigFileParser->getDocument();
	}
	catch(const XMLException& e)
	{
		string emsg("Error parsing file! :\n"+XMLString2string(e.getMessage()));
		throw ParadocksConfigError(emsg);
	}
	
	// get informations from the tree
	DOMElement* root=xmlDoc->getDocumentElement();
	if(!root)
	{
		throw ParadocksConfigError( "Error parsing file! :\nempty XML document" );
	}
	if (XMLString2string(root->getTagName())!="paradocks") 
	{
		throw ParadocksConfigError("Configuration file is not a ParaDocks config file!");
	}
	
	// get the optimizer node
	{
	XMLCh* xmloptimizer=XMLString::transcode("optimizer");
	DOMNodeList* optimizerlist=root->getElementsByTagName(xmloptimizer);
	XMLString::release(&xmloptimizer);
	const  XMLSize_t optcount=optimizerlist->getLength();
	if (optcount!=1) throw ParadocksConfigError("<optimizer> line not found in config file!");
	// get the optimizer type from Attribute type
	DOMElement* optimizerele=dynamic_cast<xercesc::DOMElement*>(optimizerlist->item(0));
	XMLCh* xmltype=XMLString::transcode("type");
	_optimizertype=XMLString2string(optimizerele->getAttribute(xmltype));
	XMLString::release(&xmltype);
	// get the par
	XMLCh* xmlpar=XMLString::transcode("par");
	DOMNodeList* optimizerparlist=optimizerele->getElementsByTagName(xmlpar);
	XMLString::release(&xmlpar);
	const  XMLSize_t optparcount=optimizerparlist->getLength();
	for(XMLSize_t i = 0; i < optparcount; ++i)
	{
		DOMElement* j=dynamic_cast<xercesc::DOMElement*>(optimizerparlist->item(i));
		XMLCh* xmlval=XMLString::transcode("val");
		_optimizerpar.push_back(XMLString2string(j->getAttribute(xmlval)));
		XMLString::release(&xmlval);
	}
	}
	
	
	// get the fitness node
	{
	XMLCh* xmlfitness=XMLString::transcode("fitness");
	DOMNodeList* fitnesslist=root->getElementsByTagName(xmlfitness);
	XMLString::release(&xmlfitness);
	const  XMLSize_t fitcount=fitnesslist->getLength();
	if (fitcount!=1) throw ParadocksConfigError("<fitness> line not found in config file!");
	// get the fitness type from Attribute type
	DOMElement* fitnessele=dynamic_cast<xercesc::DOMElement*>(fitnesslist->item(0));
	XMLCh* xmltype=XMLString::transcode("type");
	_fitnesstype=XMLString2string(fitnessele->getAttribute(xmltype));
	XMLString::release(&xmltype);
	// get the par
	XMLCh* xmlpar=XMLString::transcode("par");
	DOMNodeList* fitnessparlist=fitnessele->getElementsByTagName(xmlpar);
	XMLString::release(&xmlpar);
	const  XMLSize_t fitparcout=fitnessparlist->getLength();
	for(XMLSize_t i=0; i<fitparcout; ++i)
	{
		DOMElement* j=dynamic_cast<xercesc::DOMElement*>(fitnessparlist->item(i));
		XMLCh* xmlval=XMLString::transcode("val");
		_fitnesspar.push_back(XMLString2string(j->getAttribute(xmlval)));
		XMLString::release(&xmlval);
	}
	}
	
	// get the input node
	{
	XMLCh* xmlinput=XMLString::transcode("input");
	DOMNodeList* inputlist=root->getElementsByTagName(xmlinput);
	XMLString::release(&xmlinput);
	const  XMLSize_t inputcount=inputlist->getLength();
	if (inputcount!=1) throw ParadocksConfigError("<input> line not found in config file!");
	DOMElement* inputele=dynamic_cast<xercesc::DOMElement*>(inputlist->item(0));
	// get the protein node
	XMLCh* xmlprotein=XMLString::transcode("protein");
	DOMNodeList* proteinlist=inputele->getElementsByTagName(xmlprotein);
	XMLString::release(&xmlprotein);
	const  XMLSize_t proteincount=proteinlist->getLength();
	if (proteincount!=1) throw ParadocksConfigError("<protein> line not found in config file!");
	DOMElement* proteinele=dynamic_cast<xercesc::DOMElement*>(proteinlist->item(0));
	XMLString::release(&xmlprotein);
	XMLCh* xmlfile=XMLString::transcode("file");
	_protein=XMLString2string(proteinele->getAttribute(xmlfile));
	XMLCh* xmlx=XMLString::transcode("x");
	XMLCh* xmly=XMLString::transcode("y");
	XMLCh* xmlz=XMLString::transcode("z");
	XMLCh* xmlrad=XMLString::transcode("rad");
	_x=XMLString2string(proteinele->getAttribute(xmlx));
	_y=XMLString2string(proteinele->getAttribute(xmly));
	_z=XMLString2string(proteinele->getAttribute(xmlz));
	_rad=XMLString2string(proteinele->getAttribute(xmlrad));
	XMLString::release(&xmlx);
	XMLString::release(&xmly);
	XMLString::release(&xmlz);
	XMLString::release(&xmlrad);
	// get the ligand nodes
	XMLCh* xmlligand=XMLString::transcode("ligand");
	DOMNodeList* ligandlist=inputele->getElementsByTagName(xmlligand);
	XMLString::release(&xmlligand);
	const  XMLSize_t ligandcount=ligandlist->getLength();
	for(XMLSize_t i=0; i<ligandcount; ++i)
	{
		DOMElement* j=dynamic_cast<xercesc::DOMElement*>(ligandlist->item(i));
		_ligand.push_back(XMLString2string(j->getAttribute(xmlfile)));
		XMLCh* xmlidx=XMLString::transcode("idx");
		_idx.push_back(XMLString2string(j->getAttribute(xmlidx)));
		XMLString::release(&xmlidx);
	}
	XMLString::release(&xmlfile);
	}
	
	// get other configuration information
	{
	XMLCh* xmlconfiguration=XMLString::transcode("configuration");
	DOMNodeList* conflist=root->getElementsByTagName(xmlconfiguration);
	XMLString::release(&xmlconfiguration);
	const  XMLSize_t confcount=conflist->getLength();
	if (confcount!=1) throw ParadocksConfigError("<configuration> line not found in config file!");
	DOMElement* confele=dynamic_cast<xercesc::DOMElement*>(conflist->item(0));
	// get the configuration data
	// random seed
	XMLCh* xmlrandom=XMLString::transcode("random");
	DOMNodeList* randomlist=confele->getElementsByTagName(xmlrandom);
	XMLString::release(&xmlrandom);
	const  XMLSize_t randomcount=randomlist->getLength();
	if (randomcount!=1) throw ParadocksConfigError("<random> line not found in config file!");
	DOMElement* randomele=dynamic_cast<xercesc::DOMElement*>(randomlist->item(0));
	XMLCh* xmlseed=XMLString::transcode("seed");
	_seed=XMLString2string(randomele->getAttribute(xmlseed));
	XMLString::release(&xmlseed);
	// runs
	XMLCh* xmlruns=XMLString::transcode("runs");
	DOMNodeList* runslist=confele->getElementsByTagName(xmlruns);
	XMLString::release(&xmlruns);
	const  XMLSize_t runscount=runslist->getLength();
	if (runscount!=1) throw ParadocksConfigError("<runs> line not found in config file!");
	DOMElement* runsele=dynamic_cast<xercesc::DOMElement*>(runslist->item(0));
	XMLCh* xmlval=XMLString::transcode("val");
	_runs=XMLString2string(runsele->getAttribute(xmlval));
	XMLString::release(&xmlval);
	}
	
	// get the output node
	{
	XMLCh* xmloutput=XMLString::transcode("output");
	DOMNodeList* outputlist=root->getElementsByTagName(xmloutput);
	XMLString::release(&xmloutput);
	const  XMLSize_t outputcount=outputlist->getLength();
	if (outputcount!=1) throw ParadocksConfigError("<output> line not found in config file!");
	DOMElement* outputele=dynamic_cast<xercesc::DOMElement*>(outputlist->item(0));
	// get the output data
	// output during iteration
	XMLCh* xmliteration=XMLString::transcode("iteration");
	DOMNodeList* iterationlist=outputele->getElementsByTagName(xmliteration);
	XMLString::release(&xmliteration);
	const  XMLSize_t iterationcount=iterationlist->getLength();
	if (iterationcount!=1) throw ParadocksConfigError("<iteration> line not found in config file!");
	DOMElement* iterationele=dynamic_cast<xercesc::DOMElement*>(iterationlist->item(0));
	XMLCh* xmlval=XMLString::transcode("val");
	_iterout=XMLString2string(iterationele->getAttribute(xmlval));
	// output at the end
	XMLCh* xmlend=XMLString::transcode("end");
	DOMNodeList* endlist=outputele->getElementsByTagName(xmlend);
	XMLString::release(&xmlend);
	const  XMLSize_t endcount=endlist->getLength();
	if (endcount!=1) throw ParadocksConfigError("<end> line not found in config file!");
	DOMElement* endele=dynamic_cast<xercesc::DOMElement*>(endlist->item(0));
	_endout=XMLString2string(endele->getAttribute(xmlval));
	// prefix for output files
	XMLCh* xmlprefix=XMLString::transcode("prefix");
	DOMNodeList* prefixlist=outputele->getElementsByTagName(xmlprefix);
	XMLString::release(&xmlprefix);
	const  XMLSize_t prefixcount=prefixlist->getLength();
	if (prefixcount!=1) throw ParadocksConfigError("<prefix> line not found in config file!");
	DOMElement* prefixele=dynamic_cast<xercesc::DOMElement*>(prefixlist->item(0));
	_prefixout=XMLString2string(prefixele->getAttribute(xmlval));
	XMLString::release(&xmlval);
	}
	// shutdown XML
	try
	{
		delete ConfigFileParser;
		XMLPlatformUtils::Terminate();
	}
	catch(const XMLException& e)
	{
		string emsg("Error during XML shutdown! :\n"+XMLString2string(e.getMessage()));
		throw ParadocksConfigError(emsg);
	}
}

}
