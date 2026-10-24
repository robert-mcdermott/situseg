package situseg

// The failures a segmentation meets, classified, each with what to do about it: its prior is
// missing, its input is malformed, it would need more memory or storage than the computer has, its
// region holds nothing to segment. An error says what went wrong in a sentence, and its remedy
// what to try.

import (
	"errors"
	"fmt"
	"syscall"
)

// Failure classes.
const (
	ClassPrior   = "prior"
	ClassInput   = "input"
	ClassMemory  = "memory"
	ClassStorage = "storage"
	ClassRegion  = "region"
)

// Failed is a segmentation's classified failure.
type Failed struct {
	Class  string
	Remedy string
	Err    error
}

func (f *Failed) Error() string { return f.Err.Error() }
func (f *Failed) Unwrap() error { return f.Err }

func fail(class, remedy string, format string, args ...any) error {
	return &Failed{Class: class, Remedy: remedy, Err: fmt.Errorf(format, args...)}
}

// Remedies, by class.
const (
	remedyPrior   = "The segmentation starts from the cells the instrument gave the transcripts: name the table's cell column (-cell), and the values in it that mean no cell (-no-cell)."
	remedyInput   = "Check the table's columns (-x, -y, -gene, -cell) and their units (-scale): the places are in µm."
	remedyMemory  = "Refine a region first (-region, or the roi-pilot preset), or use the low-memory preset; a computer with more memory can refine the whole tissue."
	remedyStorage = "Free space on the disk the output is written to, or write it to a larger disk (-o)."
	remedyRegion  = "Choose a region that holds cells and transcripts: its corners are in µm, in the table's coordinates."
)

// Classify classifies an error a segmentation met: its own classified failures as they are, a full
// disk as storage, anything else unclassified (nil).
func Classify(err error) *Failed {
	var f *Failed
	if errors.As(err, &f) {
		return f
	}
	if errors.Is(err, syscall.ENOSPC) {
		return &Failed{Class: ClassStorage, Remedy: remedyStorage, Err: err}
	}
	return nil
}
