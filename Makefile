.PHONY: install check benchmark model ui
install:
	python -m pip install -e '.[dev,ui]'
check:
	python -m compileall -q pcpt grm_benchmark.py app.py
	python -m pytest -q
model:
	python -m pcpt.cli model --scenario configs/baseline.json
benchmark:
	python -m pcpt.cli benchmark
ui:
	streamlit run app.py
