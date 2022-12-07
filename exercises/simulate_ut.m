function sim_ut = simulate_ut(P, n):
    ch = dtmc(P);
    sim_ut = simulate(ch, n);
end
