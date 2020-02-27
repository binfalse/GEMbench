sed 's%/%,%' affected_subsystems.csv | \grep -v ',,,' | sed 's/EMTAB37/EMTAB-37/' > affects.csv
