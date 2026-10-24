//go:build !race

package situseg

// raceDetector: the tests run with the race detector, which multiplies memory and time.
const raceDetector = false
