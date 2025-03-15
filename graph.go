package main

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
	"math"
)

type Vertex struct {
	ID       int64
	Label    string
	Sublabel string
	OutEdges map[*Vertex]*Edge
	InEdges  map[*Vertex]*Edge
}

var counter int64 = 0

func NewVertex(label, sublabel string) *Vertex {
	v := &Vertex{
		ID:       counter,
		Label:    label,
		Sublabel: sublabel,
		OutEdges: make(map[*Vertex]*Edge, 0),
		InEdges:  make(map[*Vertex]*Edge, 0),
	}
	counter++
	return v
}

func (v *Vertex) NormalizeEdges() {
	sum := 0.0
	for _, edge := range v.OutEdges {
		sum += edge.Weight
	}
	if sum > 0 {
		for _, edge := range v.OutEdges {
			edge.Weight /= sum
		}
	}
}

type Edge struct {
	From, To *Vertex
	Weight   float64
}

func NewEdge(from, to *Vertex, weight float64) *Edge {
	return &Edge{
		From:   from,
		To:     to,
		Weight: weight,
	}
}

type Graph struct {
	Vertices []*Vertex
	Labels   map[string][]*Vertex
	Edges    map[*Vertex]map[*Vertex]*Edge
}

func NewGraph(labels map[string]int64) *Graph {
	g := &Graph{
		Vertices: make([]*Vertex, 0),
		Edges:    make(map[*Vertex]map[*Vertex]*Edge, 0),
		Labels:   make(map[string][]*Vertex, 0),
	}
	for label, num := range labels {
		g.Labels[label] = make([]*Vertex, 0)
		for i := int64(0); i < num; i++ {
			v := NewVertex(label, "")
			g.Vertices = append(g.Vertices, v)
			g.Labels[label] = append(g.Labels[label], v)
			g.Edges[v] = make(map[*Vertex]*Edge, 0)
		}
	}

	for _, u := range g.Vertices {
		for _, v := range g.Vertices {
			edge := NewEdge(u, v, 0.0)
			u.OutEdges[v] = edge
			v.InEdges[u] = edge
			g.Edges[u][v] = edge
		}
	}

	return g
}

func (g *Graph) GetEdgeWeight(from, to *Vertex) float64 {
	return g.Edges[from][to].Weight
}

func (g *Graph) SetEdgeWeight(from, to *Vertex, weight float64) {
	g.Edges[from][to].Weight = weight
}

func (g *Graph) NormalizeEdges() {
	for _, v := range g.Vertices {
		v.NormalizeEdges()
	}
}

func (g *Graph) AddVertex(label, sublabel string) *Vertex {
	v := NewVertex(label, sublabel)
	g.Vertices = append(g.Vertices, v)
	g.Edges[v] = make(map[*Vertex]*Edge)
	if _, ok := g.Labels[label]; !ok {
		g.Labels[label] = make([]*Vertex, 0)
	}
	g.Labels[label] = append(g.Labels[label], v)
	for _, u := range g.Vertices {
		edge := NewEdge(u, v, 0.0)
		u.OutEdges[v] = edge
		v.InEdges[u] = edge
		g.Edges[u][v] = edge
		edge = NewEdge(v, u, 0.0)
		v.OutEdges[u] = edge
		u.InEdges[v] = edge
		g.Edges[v][u] = edge
	}
	return v
}

func constructGraph(senders, n int64, corrupt float64, l int64, tau_ float64) {
	numCorrupt := int64(float64(n) * corrupt)
	numHonests := n - numCorrupt
	labels := make(map[string]int64)
	labels["Sender"] = senders
	labels["HonestRelay"] = numHonests
	g := NewGraph(labels)

	unCorruptNodes := make([]*Vertex, len(g.Vertices))
	for i, h := range g.Vertices {
		unCorruptNodes[i] = h
	}

	for _, h := range g.Labels["Sender"] {
		for _, j := range g.Labels["HonestRelay"] {
			g.SetEdgeWeight(h, j, 1.0)
		}
		recurse(g, h, 0, l, numCorrupt)
	}

	for _, h := range g.Labels["HonestRelay"] {
		for _, j := range g.Labels["HonestRelay"] {
			g.SetEdgeWeight(h, j, 1.0)
		}
		recurse(g, h, 0, l, numCorrupt)
	}

	g.NormalizeEdges()

	nodes := make([]*Vertex, len(g.Vertices))
	for i, h := range g.Vertices {
		nodes[i] = h
	}

	//done := 0
	//
	//for _, u := range nodes {
	//	for v, uv := range u.OutEdges {
	//		if v.Label != "Delay" && uv.Weight > 0.0 {
	//			d := g.AddVertex("Delay", "")
	//			g.SetEdgeWeight(u, d, (1.0-tau_)*uv.Weight)
	//			g.SetEdgeWeight(u, v, tau_*uv.Weight)
	//			for w, vw := range v.OutEdges {
	//				if w.Label != "Delay" && vw.Weight > 0.0 {
	//					g.SetEdgeWeight(d, w, vw.Weight)
	//				}
	//			}
	//			v.NormalizeEdges()
	//			fmt.Printf("\r             \rDone %d", done)
	//			done++
	//		}
	//	}
	//}
	//
	//fmt.Printf("\n")

	sender1 := g.Labels["Sender"][0]

	result := computeMatricesAndEvolve(g, sender1, l)

	uniform := 1.0 / float64(n)

	dif := 0.0
	for _, r := range result {
		fmt.Printf("%f\n", r)
		dif += math.Abs(r - uniform)
	}

	variationDist := dif / float64(n)

	fmt.Printf("Variation distance: %f\n", variationDist)

}

func recurse(g *Graph, root *Vertex, level int64, l, numCorrupt int64) {
	if level == l {
		return
	}
	for i := int64(0); i < numCorrupt; i++ {
		c := g.AddVertex("CorruptRelay", fmt.Sprintf("%d", i))
		g.SetEdgeWeight(root, c, 1.0)
		for _, h2 := range g.Labels["HonestRelay"] {
			g.SetEdgeWeight(c, h2, 1.0)
		}
		recurse(g, c, level+1, l, numCorrupt)
	}
}

func computeMatricesAndEvolve(g *Graph, sender *Vertex, l int64) []float64 {
	g.NormalizeEdges()

	n := len(g.Vertices) // Number of nodes

	// Create an empty matrix (n x n)
	P := mat.NewDense(n, n, nil)
	I := mat.NewVecDense(n, nil)

	for _, u := range g.Vertices {
		I.SetVec(int(u.ID), 0)
		for v, uv := range u.OutEdges {
			P.Set(int(u.ID), int(v.ID), uv.Weight)
		}
	}
	I.SetVec(int(sender.ID), 1.0)

	result := evolveProbability(P, I, int(l-1))

	resultNorm := make([]float64, 0)
	for _, h := range g.Labels["HonestRelay"] {
		resultNorm = append(resultNorm, result.AtVec(int(h.ID)))
	}
	co := make(map[string][]*Vertex)
	for _, c := range g.Labels["CorruptRelay"] {
		if _, ok := co[c.Sublabel]; !ok {
			co[c.Sublabel] = make([]*Vertex, 0)
		}
		co[c.Sublabel] = append(co[c.Sublabel], c)
	}
	for _, cs := range co {
		sum := 0.0
		for _, c := range cs {
			sum += result.AtVec(int(c.ID))
		}
		resultNorm = append(resultNorm, sum)
	}
	return resultNorm
}

// Function to evolve the probability distribution
func evolveProbability(P *mat.Dense, x0 *mat.VecDense, steps int) *mat.VecDense {
	state := mat.NewVecDense(x0.Len(), nil)
	state.CopyVec(x0) // Initialize with x0

	for t := 0; t < steps; t++ {
		nextState := mat.NewVecDense(state.Len(), nil)
		nextState.MulVec(P, state) // x_{t+1} = P * x_t

		// Print the state vector at each step
		fmt.Printf("Step %d: %v\n", t+1, nextState.RawVector().Data)

		// Update state
		state.CopyVec(nextState)
	}

	return state
}

//
//type CustomNode struct {
//	graph.Node  // Embeds gonum's Node interface
//	Label       string
//	IsCorrupted bool
//	IsSender    bool
//}
//
//var nodes = make(map[int64]*CustomNode)
//var g *simple.WeightedDirectedGraph = simple.NewWeightedDirectedGraph(0, 0) // No negative weights
//var nodeLabels = make(map[string][]*CustomNode)
//var tau float64 = 0.0
//
//var counter int64 = 0
//
//func addNode(label string, iscorrupted bool, issender bool) *CustomNode {
//	node := &CustomNode{
//		Node:        simple.Node(counter),
//		Label:       label,
//		IsCorrupted: iscorrupted,
//		IsSender:    issender,
//	}
//	counter++
//	nodes[node.ID()] = node
//	g.AddNode(node)
//	if _, ok := nodeLabels[label]; !ok {
//		nodeLabels[label] = make([]*CustomNode, 0)
//	} else {
//		nodeLabels[label] = append(nodeLabels[label], node)
//	}
//	return node
//}
//
//func getNodesByLabel(label string) []*CustomNode {
//	return nodeLabels[label]
//}
//
//func removeNode(id int64) bool {
//	if node, ok := nodes[id]; !ok {
//		pl.LogNewError("Node with id %d not found.", id)
//		return false
//	} else {
//		g.RemoveNode(node.ID())
//		return true
//	}
//}
//
//func getNode(id int64) (*CustomNode, bool) {
//	node, ok := nodes[id]
//	return node, ok
//}
//
//func addEdge(from, to int64, weight float64) {
//	if fromNode, ok := getNode(from); !ok {
//		pl.LogNewError("Node with id %d not found.", from)
//		return
//	} else if toNode, ok2 := getNode(to); !ok2 {
//		pl.LogNewError("Node with id %d not found.", to)
//		return
//	} else {
//		if we := g.WeightedEdge(from, to); we == nil {
//			g.SetWeightedEdge(g.NewWeightedEdge(fromNode, toNode, weight))
//		} else {
//			g.SetWeightedEdge(g.NewWeightedEdge(fromNode, toNode, weight+we.Weight()))
//		}
//	}
//}
//
//func getEdge(from, to int64) (float64, bool) {
//	if we := g.WeightedEdge(from, to); we != nil {
//		return we.Weight(), true
//	} else {
//		return 0, false
//	}
//}
//
//func removeEdge(from, to int64, weight float64) bool {
//	if _, ok := getNode(from); !ok {
//		pl.LogNewError("Node with id %d not found.", from)
//		return false
//	} else if _, ok2 := getNode(to); !ok2 {
//		pl.LogNewError("Node with id %d not found.", to)
//		return false
//	} else {
//		g.RemoveEdge(from, to)
//		return true
//	}
//}
//
//// Create a directed weighted graph
//func createGraph(senders, n int64, corrupt float64, l int64, tau_ float64) *simple.WeightedDirectedGraph {
//
//	tau = tau_
//
//	for i := int64(0); i < senders; i++ {
//		addNode("Sender", false, true)
//	}
//
//	numCorrupt := int64(float64(n) * corrupt)
//	numHonests := n - numCorrupt
//
//	for i := int64(0); i < numHonests; i++ {
//		addNode("HonestRelay", false, false)
//	}
//
//	nonCorrupt := make([]*CustomNode, senders+numHonests)
//	for i, h := range getNodesByLabel("HonestRelay") {
//		nonCorrupt[i] = h
//	}
//	for i, h := range getNodesByLabel("Sender") {
//		nonCorrupt[i+int(numHonests)] = h
//	}
//
//	for _, h := range nonCorrupt {
//		for _, j := range getNodesByLabel("HonestRelay") {
//			addEdge(h.ID(), j.ID(), 1.0)
//		}
//		recurse(h.ID(), 0, l, numCorrupt)
//	}
//
//	for edge, _ := range g.WeightedEdges() {
//
//	}
//
//
//
//	return g
//}
//
//func recurse(root int64, level int64, l, numCorrupt int64) {
//	if level == l {
//		return
//	}
//	for i := int64(0); i < numCorrupt; i++ {
//		c := addNode(fmt.Sprintf("CorruptRelay_%d", i), true, false)
//		addEdge(root, c.ID(), 1.0)
//		for _, h2 := range getNodesByLabel("HonestRelay") {
//			addEdge(c.ID(), h2.ID(), 1.0)
//		}
//		recurse(c.ID(), level+1, l, numCorrupt)
//	}
//}
