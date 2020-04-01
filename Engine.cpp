#include "Engine.h"
#include <string>
#include <iostream>
#include <algorithm>

void Engine::init_system(const char* fname)
{
    std::ifstream fparticle{ fname };
    while (fparticle.peek() == '#') {
        std::string type;
        fparticle >> type;
        if (type == "#gravity:") {
            fparticle >> G.x() >> G.y() >> G.phi();
            fparticle.ignore(100, '\n');
            std::cout << "gravity: " << G << std::endl;
        }
        else if (type == "#Time:") {
            fparticle >> Time;
            fparticle.ignore(100, '\n');
            std::cout << "Time: " << Time << std::endl;
        }
        else if (type == "#nstep:") {
            fparticle >> nstep;
            fparticle.ignore(100, '\n');
            std::cout << "nstep: " << nstep << std::endl;
        }
        else if (type == "#timestep:") {
            fparticle >> timestep;
            fparticle.ignore(100, '\n');
            std::cout << "timestep: " << timestep << std::endl;
        }
        else if (type == "#nprint:") {
            fparticle >> nprint;
            fparticle.ignore(100, '\n');
            std::cout << "nprint: " << nprint << std::endl;
        }
        else if (type == "#nenergy:") {
            fparticle >> nenergy;
            fparticle.ignore(100, '\n');
            std::cout << "nenergy: " << nenergy << std::endl;
        }
        else if (type == "#lx:") {
            fparticle >> lx;
            fparticle.ignore(100, '\n');
            std::cout << "lx: " << lx << std::endl;
        }
        else if (type == "#ly:") {
            fparticle >> ly;
            fparticle.ignore(100, '\n');
            std::cout << "ly: " << ly << std::endl;
        }
        else if (type == "#x_0:") {
            fparticle >> x_0;
            fparticle.ignore(100, '\n');
            std::cout << "x_0: " << x_0 << std::endl;
        }
        else if (type == "#y_0:") {
            fparticle >> y_0;
            fparticle.ignore(100, '\n');
            std::cout << "y_0: " << y_0 << std::endl;
        }
        else {
            std::cerr << "init: unknown global property: " << type << std::endl;
            abort();
        }
    }
    while (fparticle) {
        Sphere pp;
        fparticle >> pp;
        if (fparticle) {
            particle.push_back(pp);
        }
    }
    no_of_particles = particle.size();
    std::cout << no_of_particles << " particles read" << std::endl;
}

void Engine::step()
{
    integrate();
    check_dump();
}

void Engine::phase_plot(std::ostream& os)
{
    os << "#NewFrame\n";
    os << "#no_of_particles: " << no_of_particles << std::endl;
    os << "#compressed: no\n";
    os << "#type: SphereXYPhiVxVyOmegaRMFixed25\n";
    os << "#gravity: " << G.x() << " " << G.y() << " " << G.phi() << std::endl;
    os << "#Time: " << Time << std::endl;
    os << "#timestep: " << timestep << std::endl;
    os << "#EndOfHeader\n";
    for (unsigned int i = 0; i < particle.size(); i++) {
        os << particle[i];
    }
    os << std::flush;
}

double Engine::total_kinetic_energy()
{
    double sum = 0;
    for (unsigned int i{ 0 }; i < particle.size(); i++) {
        if (particle[i].ptype() == 0) {
            sum += particle[i].kinetic_energy();
        }
    }
    return sum;
}

std::vector<Vector> Engine::get_particle_positions()
{
    std::vector<Vector> result;
    for (auto p : particle) {
        Vector vec{ p.x(), p.y(), p.r() };
        result.push_back(vec);
    }
    return result;
}

void Engine::make_forces()
{
    for (unsigned int i{ 0 }; i < no_of_particles - 1; i++) {
        for (unsigned int k{ i + 1 }; k < no_of_particles; k++) {
            if ((particle[i].ptype() == 0) || (particle[k].ptype() == 0)) {
                force(particle[i], particle[k], lx, ly);
            }
        }
    }
}

void Engine::check_dump()
{
    if (save != save_interval) { 
        save++; 
    }
    else {
        save = 0;
        dump();
    }
}


void Engine::dump()
{
    std::fprintf(f1, "ITEM: TIMESTEP\n%d\n", int(Time / timestep));
    std::fprintf(f1, "ITEM: BOX BOUNDS pp pp pp\n%.4f %.4f\n%.4f %.4f\n0 1\n", x_0, x_0 + lx, y_0, y_0 + ly);
    std::fprintf(f1, "ITEM: NUMBER OF ATOMS\n%d\n", no_of_particles);
    std::fprintf(f1, "ITEM: ATOMS x y z radius type\n");
    for (Sphere& p : particle) {
        std::fprintf(f1, "%.3f %.3f %.3f %.3f %d\n", p.x(), p.y(), 0.0, p.r(), p.ptype());
    }
}

void Engine::integrate()
{
    std::for_each(particle.begin(), particle.end(),
        [&](Sphere& p) {
            if (p.ptype() == 0) {
                p.set_force_to_zero();
                p.predict(timestep);
            }
            else {
                p.boundary_conditions(timestep, Time);
            }
        });

    make_forces();

    std::for_each(particle.begin(), particle.end(),
        [&](Sphere& p) {if (p.ptype() == 0) p.correct(timestep, G); });

    std::for_each(particle.begin(), particle.end(),
        [&](Sphere& p) {p.periodic_bc(x_0, y_0, lx, ly); });

    Time += timestep;
}


