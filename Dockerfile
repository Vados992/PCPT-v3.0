FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1
WORKDIR /app
COPY pyproject.toml README.md LICENSE grm_benchmark.py app.py ./
COPY pcpt/ ./pcpt/
RUN pip install --no-cache-dir '.[ui]'
COPY configs/ ./configs/
COPY data/templates/ ./data/templates/
RUN useradd -u 10001 -m pcpt && mkdir -p /app/outputs && chown -R pcpt:pcpt /app/outputs
USER pcpt
EXPOSE 8501
CMD ["streamlit", "run", "app.py", "--server.address=0.0.0.0", "--server.port=8501", "--server.headless=true"]
