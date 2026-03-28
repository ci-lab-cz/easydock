FROM dptechnology/unimol:latest-pytorch1.11.0-cuda11.3

COPY unipka.py /unipka/unipka.py
COPY simple_smarts_pattern.tsv /unipka/simple_smarts_pattern.tsv
COPY smarts_pattern.tsv /unipka/smarts_pattern.tsv

RUN wget -q -O /unipka/t_dwar_v_novartis_a_b.pt \
    https://zenodo.org/records/19288649/files/t_dwar_v_novartis_a_b.pt

ENTRYPOINT ["/bin/sh", "-c", "\
case \"$1\" in \
  protonate) shift; exec python /unipka/unipka.py \"$@\" ;; \
  help) echo \"Use docker run <image> protonate <args>\" ;; \
  *) exec \"$@\" ;; \
esac", "--"]