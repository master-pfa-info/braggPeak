package main

import (
	"fmt"
	"math"

	"go-hep.org/x/hep/hplot"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

const (
	re = 2.818e-13 // electron classical radius (cm)
	me = 0.511     // electron mass (MeV)
	mp = 938.2720  // proton mass (MeV)
)

type Particle struct {
	M  float64 // mass (MeV)
	T  float64 // kinetic energy (MeV)
	T0 float64 // initial kinetic energy (MeV)
	Z  float64 // charge (in units of e)
}

func NewParticle(m, t0, z float64) *Particle {
	return &Particle{M: m, T: t0, T0: t0, Z: z}
}

func (p *Particle) BetaGamma() (float64, float64) {
	beta := math.Sqrt(1 - p.M*p.M/(p.M+p.T)/(p.M+p.T))
	gamma := 1 / math.Sqrt(1-beta*beta)
	return beta, gamma
}

func (p *Particle) Momentum() float64 {
	return math.Sqrt(p.T * (p.T + 2*p.M))
}

type Material struct {
	N float64 // number of electrons per unit volume (in cm^{-3})
	I float64 // mean ionization potential (MeV)
}

func NewMaterial(n, i float64) *Material {
	return &Material{N: n, I: i}
}

func dEdx(p *Particle, m *Material) float64 {
	beta, gamma := p.BetaGamma()
	betagamma := beta * gamma
	constant := 4 * math.Pi * me * re * re
	F := math.Log(2*me*betagamma*betagamma) - beta*beta
	//dEdx := 4 * math.Pi * me * re * re * m.Z * m.Rho / m.A * Na * p.Z * p.Z / (beta * beta) * (math.Log(2*me*betagamma*betagamma/m.I) - beta*beta)
	dEdx := constant * m.N * p.Z * p.Z / (beta * beta) * (F - math.Log(m.I))
	fmt.Println("Debug:", constant, F, p.T, p.Momentum(), beta, gamma, dEdx)
	return dEdx
}

type BraggPeak struct {
	Ts  []float64
	dEs []float64
	Xs  []float64
}

func NewBraggPeak(p *Particle, m *Material) *BraggPeak {
	bp := &BraggPeak{}
	dx := 0.01e-1 // in cm

	var Ts []float64
	var dEs []float64
	var Xs []float64
	n := 0
	for p.T >= 0 {
		n++

		dEdx := dEdx(p, m)
		dE := dEdx * dx
		p.T -= dE
		dEs = append(dEs, dE)
		Ts = append(Ts, p.T)
		Xprev := 0.
		if len(Xs) > 0 {
			Xprev = Xs[len(Xs)-1]
		}
		Xs = append(Xs, Xprev+dx)
	}

	bp.Ts = Ts
	bp.dEs = dEs
	bp.Xs = Xs
	return bp
}

func (b *BraggPeak) Len() int {
	return len(b.Ts)
}

func (b *BraggPeak) XY(i int) (float64, float64) {
	return b.Xs[i], b.dEs[i]
}

func (b *BraggPeak) Range() float64 {
	return b.Xs[len(b.Xs)-1]
}

func main() {
	// 	materialH := NewMaterial(1*1/9., 1, 1, 20*10e-6)
	// 	materialO := NewMaterial(1*8/9., 16, 8, 8*12*10e-6)
	material := NewMaterial(3.34e23, 78e-6)

	/*
		Ts := make([]float64, 10000)
		dEdxs := make([]float64, len(Ts))
		for i := range Ts {
			Ts[i] = 1 + float64(i)*(1e5-1)/float64(len(Ts))
			proton := NewParticle(mp, Ts[i], 1)
			dEdxs[i] = dEdx(proton, material)
		}

		pts := make(plotter.XYs, len(Ts))
		for i := range Ts {
			pts[i].X = Ts[i]
			pts[i].Y = dEdxs[i]
		}
		p, err := plot.New()
		if err != nil {
			panic(err)
		}

		p.Title.Text = ""
		p.X.Label.Text = "T"
		p.Y.Label.Text = "-dE/dx"
		//p.X.Tick.Marker = &hplot.FreqTicks{N: 10, Freq: 1}
		//p.X.Scale = plot.LogScale{}
		//p.Y.Scale = plot.LogScale{}
		p.Add(hplot.NewGrid())
		err = plotutil.AddScatters(p, pts)
		if err != nil {
			panic(err)
		}
		// Save the plot to a PNG file.
		if err := p.Save(6*vg.Inch, 3*vg.Inch, "StoppingPower.png"); err != nil {
			panic(err)
		}
	*/

	// The value of I is taken from http://www.sciencedirect.com/science/article/pii/S135044870700409X

	proton := NewParticle(mp, 60, 1)
	fmt.Println(dEdx(proton, material))

	bp65 := NewBraggPeak(NewParticle(mp, 65, 1), material)
	bp100 := NewBraggPeak(NewParticle(mp, 100, 1), material)
	bp150 := NewBraggPeak(NewParticle(mp, 150, 1), material)
	bp200 := NewBraggPeak(NewParticle(mp, 200, 1), material)

	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Title.Text = ""
	p.X.Label.Text = "X"
	p.Y.Label.Text = "T"
	//p.X.Tick.Marker = &hplot.FreqTicks{N: 10, Freq: 1}
	//p.X.Scale = plot.LogScale{}
	//p.Y.Scale = plot.LogScale{}
	p.Add(hplot.NewGrid())
	err = plotutil.AddScatters(p, bp65, bp100, bp150, bp200)
	if err != nil {
		panic(err)
	}
	// Save the plot to a PNG file.
	if err := p.Save(6*vg.Inch, 3*vg.Inch, "dEvsX.png"); err != nil {
		panic(err)
	}

	fmt.Println(bp65.Range(), bp100.Range(), bp150.Range(), bp200.Range())

}
