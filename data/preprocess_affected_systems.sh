sed 's%/%,%' affected_subsystems.csv | \grep -v ',,,' > affects.csv
