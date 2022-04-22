FROM python:3.10.4-slim
LABEL maintainer="The NuRadioReco Authors <physics-astro-nuradiomcdev@lists.uu.se>"

RUN apt-get update
RUN apt-get upgrade -y

# copy NuRadioMC
COPY . /usr/local/lib/NuRadioMC
# Install core dependencies
RUN python /usr/local/lib/NuRadioMC/install_dev.py --install --no-interactive

# Install optional dependencies #TODO: add to pyproject.toml
RUN pip install dnspython gunicorn DateTime pandas

#Uninstall and reinstall werkzeug bug
#RUN pip uninstall Werkzeug
RUN pip install Werkzeug==2.0.0

# Add NuRadioMC to PYTHONPATH
ENV PYTHONPATH="/usr/local/lib/NuRadioMC:$PYTHONPATH"

RUN useradd nuradio

USER nuradio
EXPOSE 8050
WORKDIR /usr/local/lib/NuRadioMC/NuRadioReco/detector/webinterface
CMD python index.py