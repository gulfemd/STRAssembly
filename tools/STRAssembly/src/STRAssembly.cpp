//! [snippet1]

#include <STRAssembly.hpp>
#include <gatb/debruijn/impl/Simplifications.hpp>

using namespace std;

/********************************************************************************/

// We define some constant strings for names of command line parameters
static const char* STR_FOO = "-foo";
static const char* STR_SAM_FILE_PATH = "-sam";
static const char* STR_OUTPUT_DIR_PATH = "-out";

string processBamRecord (string line)
{
    stringstream ss (line);
    string qname, flag, rname, pos, mapQ, cigar, rnext, pnext, tlen, seq, qual;
    ss >> qname >> flag >> rname >> pos >> mapQ >> cigar >> rnext >> pnext >> tlen >> seq >> qual;
    return seq;
}

bool discardSequence (string seq)
{
    size_t found = seq.find('.');
    if (found != string::npos)
        return false;

    found = seq.find('N');
    if (found != string::npos)
        return false;

    return true;
}

string baseName(const string &path)
{
    return path.substr(path.find_last_of("/\\") + 1);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void STRAssembly::assemble (GraphType &graph)
{
    string samFile = baseName(getInput()->getStr(STR_SAM_FILE_PATH));
    string outputPath = getInput()->getStr(STR_OUTPUT_DIR_PATH);

    if (outputPath.back() != '/')
        outputPath += '/';

    outputPath += samFile.substr(0, samFile.size() - 4) + ".fa";

    IBank* outputBank = new BankFasta (outputPath);

    GraphIterator<NodeGU> it = graph.iterator ();
    cout << "Number of nodes: " << it.size() << endl;

    bool verbose = false;
    Simplifications<GraphType,NodeGU,EdgeGU> simplifier(graph, 1, verbose);

    unsigned long nbTipsRemoved = simplifier.removeTips();
    cout << "nbTipsRemoved: " << nbTipsRemoved << endl;

    // unsigned long nbBulgesRemoved = simplifier.removeBulges();
    // cout << "nbBulgesRemoved: " << nbBulgesRemoved << endl;

    // unsigned long nbECsRemoved = simplifier.removeErroneousConnections();
    // cout << "nbECsRemoved: " << nbECsRemoved << endl;

    GraphIterator<NodeGU> it2 = graph.iterator ();
    cout << "Number of nodes: " << it2.size() << endl;

    ProgressGraphIteratorTemplate<NodeGU,ProgressTimerAndSystem> itNode (graph.GraphType::iterator(), "str : assembly");

    for (itNode.first(); !itNode.isDone(); itNode.next()) {

        NodeGU node = itNode.item();

        if (graph.unitigIsMarked(node)) {
            continue;
        }

        if (graph.isNodeDeleted(node)) {
            continue;
        }

        cout << endl << "-------------------------- " << graph.toString (node) << " -------------------------" << endl;

        assembleFrom(node, graph, outputBank);
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void STRAssembly::assembleFrom (NodeGU startingNode, GraphType& graph, IBank *outputBank)
{
    bool keepIsolatedTigs = false;
    bool isolatedLeft, isolatedRight;

    float coverage = 0;
    unsigned int isolatedCutoff = std::max(2*(unsigned int)graph.getKmerSize(), (unsigned int)150);

    string sequence = graph.simplePathBothDirections(startingNode, isolatedLeft, isolatedRight, true, coverage);

    Sequence seq (Data::ASCII);
    seq.getData().setRef ((char*)sequence.c_str(), sequence.size());
    /** We set the sequence comment. */
    stringstream ss1;
    // spades-like header (compatible with bandage)
    ss1 << "NODE_"<< nbContigs + 1 << "_length_" << sequence.size() << "_cov_" << fixed << std::setprecision(3) << coverage << "_ID_" << nbContigs;
    seq._comment = ss1.str();

    unsigned int lenTotal = sequence.size();

    if (lenTotal > isolatedCutoff || (lenTotal <= isolatedCutoff && (!(isolatedLeft && isolatedRight))) || keepIsolatedTigs) {
        outputBank->insert (seq);
        nbContigs += 1;
        totalNt += lenTotal;

        if (lenTotal > maxContigLen) {
            maxContigLen = lenTotal;
        }
    }
    else {
        nbSmallContigs++;
    }

    return;
}


/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
STRAssembly::STRAssembly ()  : Tool ("STRAssembly")
{
    // We add some custom arguments for command line interface
    getParser()->push_front (new OptionOneParam (STR_SAM_FILE_PATH, "Path for SAM file containing reads.", true));
    getParser()->push_front (new OptionOneParam (STR_OUTPUT_DIR_PATH, "Path for output directory.", true));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void STRAssembly::execute ()
{
    // We can do here anything we want.
    // For further information about the Tool class, please have a look
    // on the ToyTool snippet  (http://gatb-core.gforge.inria.fr/snippets_tools.html)
    string samPath = getInput()->getStr(STR_SAM_FILE_PATH);
    ifstream in(samPath);
    string line;
    vector<string> sequencesData;

    while (getline (in, line)) {
        if (line[0] == '@')
            continue;
        string sequence = processBamRecord (line);
        if (discardSequence (sequence))
            sequencesData.push_back (sequence);
    }

    cout << "Vector size: " << sequencesData.size() << endl;

    size_t kmerSize = 63;

    GraphType graph = GraphType::create(new BankStrings(sequencesData), "-kmer-size %d  -abundance-min 1  -verbose 0", kmerSize);

    assemble(graph);
}
