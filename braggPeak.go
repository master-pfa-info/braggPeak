package main

import (
	"math"

	"go-hep.org/x/hep/hplot"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

const (
	re = 2.818e-13 // electron classical radius (cm)
	me = 0.511     // electron mass (MeV)
	Na = 6.022e23  // Avogadro´s constant
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
	Rho float64 // volumic mass (g/cm^3)
	A   float64 // mass number
	Z   float64 // atomic number
	I   float64 // mean ionization potential (MeV)
}

func NewMaterial(rho, a, z, i float64) *Material {
	return &Material{Rho: rho, A: a, Z: z, I: i}
}

func dEdx(p *Particle, m *Material) float64 {
	beta, gamma := p.BetaGamma()
	betagamma := beta * gamma
	dEdx := 4 * math.Pi * me * re * re * m.Z * m.Rho / m.A * Na * p.Z * p.Z / (beta * beta) * (math.Log(2*me*betagamma*betagamma/m.I) - 2*beta*beta)
	//fmt.Println("Debug:", p.T, p.Momentum(), beta, gamma, dEdx)
	return dEdx
}

func main() {
	/*
		material := NewMaterial(1, 16, 8, 8*10e-6)
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

	material := NewMaterial(1, 4, 2, 6*13e-6)
	proton := NewParticle(mp, 140, 1)
	dx := 0.1e-1 // in cm

	var Ts []float64
	var dEs []float64
	var Xs []float64
	n := 0
	for proton.T >= 0 {
		n++
		dEdx := dEdx(proton, material)
		dE := dEdx * dx
		proton.T -= dE
		dEs = append(dEs, dE)
		Ts = append(Ts, proton.T)
		Xs = append(Xs, float64(n)*dx)
	}

	pts := make(plotter.XYs, len(Ts))
	for i := range Ts {
		pts[i].X = Xs[i]
		pts[i].Y = dEs[i]
	}
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
	err = plotutil.AddScatters(p, pts)
	if err != nil {
		panic(err)
	}
	// Save the plot to a PNG file.
	if err := p.Save(6*vg.Inch, 3*vg.Inch, "TvsX.png"); err != nil {
		panic(err)
	}

}
