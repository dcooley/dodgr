#include "astar.h"
#include "dgraph.h"
#include "heaps/heap.h"

#include <algorithm> // std::fill

Astar::Astar(unsigned int n, const HeapDesc& heapD, std::shared_ptr<const DGraph> g) 
{
    m_heap = heapD.newInstance(n);
    m_s = new bool[n];
    m_f = new bool[n];
    init(g);
}

/* - Destructor - */
Astar::~Astar() {
    delete [] m_s;
    delete [] m_f;
    delete m_heap;
}

/* - init() -
 * Initialise the algorithm for use with the graph pointed to by g.
 */
void Astar::init(std::shared_ptr<const DGraph> g) {
    m_graph = g;
}

/* - run() -
 * Run the algorithm, computing single-source from the starting vertex v0.
 * This assumes that the array d has been initialised with d[v] = INFINITE_DIST
 * for all vertices v != v0.
 */

void Astar::run (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<int>& prev,
        unsigned int v0)
{
    /* indexes, counters, pointers */
    const DGraphEdge *edge;

    /*** initialisation ***/

    /* optimise access to the data structures allocated for the algorithm */
    const unsigned int n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    /* initialise all vertices as unexplored */
    std::fill (m_s, m_s + n, false);
    std::fill (m_f, m_f + n, false);

    /* heap is used for heuristic dists; final dists are here: */
    std::vector <double> d_final (n);

    /* place v0 into the frontier set with a distance of zero */
    w [v0] = 0.0;
    d [v0] = 0.0;
    prev [v0] = -1;
    m_heap->insert(v0, 0.0);
    m_f [v0] = true;

    /* repeatedly update distances from the minimum remaining trigger vertex */
    while (m_heap->nItems() > 0) {
        /* delete the vertex in frontier that has minimum distance */
        unsigned int v = m_heap->deleteMin();

        /* the selected vertex moves from the frontier to the solution set */
        m_s [v] = true;
        m_f [v] = false;

        /* explore the OUT set of v */
        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!m_s [et]) {
                double wt = w [v] + edge->wt;
                if (wt < w [et]) {
                    d [et] = d [v] + edge->dist;
                    w [et] = wt;
                    prev [et] = static_cast <int> (v);
                    if (m_f [et]) {
                      m_heap->decreaseKey(et, wt);
                    }
                    else {
                      m_heap->insert (et, wt);
                      m_f [et] = true;
                    }
                }
            }

            edge = edge->nextOut;
        } /* while */
    } /* while */
}
