package main

import (
	pl "github.com/HannahMarsh/PrettyLogger"
)

func main() {
	// Create graph

	pl.SetUpLogrusAndSlog("debug")

	constructGraph(200, 100, 0.5, 3, 0.8)
}
