FROM python:3.9-slim

# Set working directory
WORKDIR /app

# Copy requirements first (for better caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the code
COPY . .

# Install the package
RUN pip install -e .

# Command to run the script
ENTRYPOINT ["python", "src/main.py"]
