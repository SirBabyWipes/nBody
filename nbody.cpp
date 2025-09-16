#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include <chrono>


const double G = 6.674e-11;
const double softing = 0.0000001;

class Particle {
public:
    double mass;
    // (x, y, z)
    std::vector<double> position;
    std::vector<double> velocity;
    std::vector<double> force;

    Particle(double m, std::vector<double> p, std::vector<double> v) {
        mass = m;
        position = p;
        velocity = v;
        force = {0, 0, 0};
    }
};

// random initialization
std::vector<Particle> iRandom(int n){
    std::vector<Particle> particles;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> rMass(1e10, 1e15);
    std::uniform_real_distribution<> rPosition(-1000, 1000);
    std::uniform_real_distribution<> rVelocity(-1, 1);

    for(int i = 0; i < n; i++){
        double m = rMass(gen);
        std::vector<double> p = {rPosition(gen), rPosition(gen), rPosition(gen)};
        std::vector<double> v = {rVelocity(gen), rVelocity(gen), rVelocity(gen)};
        particles.push_back(Particle(m, p, v));
    }
    return particles;
}

// predefined initialization
std::vector<Particle> iPredefined(){
    std::vector<Particle> particles;

    Particle sun(100, {0, 0, 0}, {0, 0, 0});
    particles.push_back(sun);

    Particle earth(10, {100, 0, 0}, {100, 0, 0});
    particles.push_back(earth);

    Particle moon(1, {105, 0, 0}, {100, 0, 0});
    particles.push_back(moon);
    
    return particles;
}

// file initialization
std::vector<Particle> iFile(const std::string &filename){
    std::vector<Particle> particles;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return particles;
    }

    int n;  
    file >> n;

    for (int i = 0; i < n; i++) {
        double m, x, y, z, vx, vy, vz, fx, fy, fz;
        file >> m >> x >> y >> z >> vx >> vy >> vz >> fx >> fy >> fz;
        if (!file) break;

        Particle p(m, {x, y, z}, {vx, vy, vz});
        p.force = {fx, fy, fz};

        particles.push_back(p);
    }

    return particles;
}

void forces(std::vector<Particle> &particles){
    // set force to 0
    for( Particle &p: particles){
        p.force = {0, 0, 0};
    }
    // loop over all particles per particle
    for(int i = 0; i < particles.size(); i++){
        for (int j = i + 1; j < particles.size(); j++){
            // calculate distance between particles
            double dx = particles[j].position[0] - particles[i].position[0];
            double dy = particles[j].position[1] - particles[i].position[1];
            double dz = particles[j].position[2] - particles[i].position[2];
            double distance = sqrt((dx*dx)+(dy*dy)+(dz*dz) + softing*softing);

            // calculate gravitational force
            double force = G * ((particles[i].mass * particles[j].mass) / (distance * distance));

            // apply force to particles
            particles[i].force[0] += force * (dx / distance);
            particles[i].force[1] += force * (dy / distance);
            particles[i].force[2] += force * (dz / distance);

            particles[j].force[0] -= force * (dx / distance);
            particles[j].force[1] -= force * (dy / distance);
            particles[j].force[2] -= force * (dz / distance);
        }
    }
}

void updateVelocity(std::vector<Particle> &particles, double dt){
    for (Particle &p : particles) {
        p.velocity[0] += (p.force[0] / p.mass) * dt;
        p.velocity[1] += (p.force[1] / p.mass) * dt;
        p.velocity[2] += (p.force[2] / p.mass) * dt;
    }
}

void updatePosition(std::vector<Particle> &particles, double dt){
    for (auto &p : particles) {
        p.position[0] += p.velocity[0] * dt;
        p.position[1] += p.velocity[1] * dt;
        p.position[2] += p.velocity[2] * dt;
    }
}

int main(int argc, char* argv[]) {
    double dt;
    int steps;
    int dumpInterval;
    std::vector<Particle> particles;

    if (argc != 5) {
        dt = 200;
        steps = 5000000;
        dumpInterval = 100;
        particles = iFile("solar.tsv");
    } else {
        int n = std::stoi(argv[1]);
        dt = std::stod(argv[2]);
        steps = std::stoi(argv[3]);
        dumpInterval = std::stoi(argv[4]);
        particles = iRandom(n);
    }
    

    // predefined initialization
    // particles = iPredefined();

    /* check step 1
    std::cout << "Loaded " << particles.size() << " particles:\n";
    for (size_t i = 0; i < particles.size(); i++) {
        const Particle &p = particles[i];
        std::cout << "Particle " << i << ":\n";
        std::cout << "  Mass = " << p.mass << "\n";
        std::cout << "  Position = (" 
                  << p.position[0] << ", " 
                  << p.position[1] << ", " 
                  << p.position[2] << ")\n";
        std::cout << "  Velocity = (" 
                  << p.velocity[0] << ", " 
                  << p.velocity[1] << ", " 
                  << p.velocity[2] << ")\n";
        std::cout << "  Force = (" 
                  << p.force[0] << ", " 
                  << p.force[1] << ", " 
                  << p.force[2] << ")\n";
    }
    */

    /* check step 2
    std::cout << "Forces after calculation:\n";
    for (size_t i = 0; i < particles.size(); i++) {
        const Particle &p = particles[i];
        std::cout << "Particle " << i 
                  << " Force = (" 
                  << p.force[0] << ", " 
                  << p.force[1] << ", " 
                  << p.force[2] << ")\n";
    }
    */

    // console out
    /*
    for (int step = 0; step < steps; step++) {
        forces(particles);
        updateVelocity(particles, dt);
        updatePosition(particles, dt);

        std::cout << "Step " << step << ":\n";
        for (size_t i = 0; i < particles.size(); i++) {
            const Particle &p = particles[i];
            std::cout << "Particle " << i
                      << " Position = (" << p.position[0] << ", "
                                         << p.position[1] << ", "
                                         << p.position[2] << ") "
                      << " Velocity = (" << p.velocity[0] << ", "
                                         << p.velocity[1] << ", "
                                         << p.velocity[2] << ")\n";
        }
    }
    */
    // file out
    std::ostringstream fname;
    fname << "sim_output_" << particles.size() << ".tsv";
    std::ofstream out(fname.str());

    // start time
    auto start = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < steps; step++) {
        forces(particles);
        updateVelocity(particles, dt);
        updatePosition(particles, dt);

        // dumpstate at intervals
        if (step % dumpInterval == 0) {
            out << particles.size();
            for (const Particle &p : particles) {
                out << "\t" << p.mass
                    << "\t" << p.position[0]
                    << "\t" << p.position[1]
                    << "\t" << p.position[2]
                    << "\t" << p.velocity[0]
                    << "\t" << p.velocity[1]
                    << "\t" << p.velocity[2]
                    << "\t" << p.force[0]
                    << "\t" << p.force[1]
                    << "\t" << p.force[2];
            }
            out << "\n";
        }
    }

    out.close();

    // end and print time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> timeElapsed = end - start;
    std::cout << "Simulation completed in " << timeElapsed.count() << " seconds.\n";

    return 0;
}
