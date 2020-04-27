FROM ubuntu:18.04

RUN apt-get update && apt-get install -y openssh-server git \
    unzip wget build-essential zlib1g-dev make iproute2 cmake \
    python python-pip openmpi-bin openmpi-common libopenmpi-dev iputils-ping

RUN git clone https://github.com/marchartung/ctsat ctsat
RUN cd ctsat && ./build_mpi.sh

ADD * ctsat/

ADD mpi-run.sh supervised-scripts/mpi-run.sh
ADD make_combined_hostfile.py supervised-scripts/make_combined_hostfile.py
RUN cd ..

RUN chmod 755 supervised-scripts/mpi-run.sh


RUN pip install supervisor awscli

RUN mkdir /var/run/sshd
RUN echo 'root:THEPASSWORDYOUCREATED' | chpasswd
RUN sed -i 's/PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

RUN echo "${USER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
ENV SSHDIR /root/.ssh
RUN mkdir -p ${SSHDIR}
RUN touch ${SSHDIR}/sshd_config
RUN ssh-keygen -t rsa -f ${SSHDIR}/ssh_host_rsa_key -N ''
RUN cp ${SSHDIR}/ssh_host_rsa_key.pub ${SSHDIR}/authorized_keys
RUN cp ${SSHDIR}/ssh_host_rsa_key ${SSHDIR}/id_rsa
RUN echo " IdentityFile ${SSHDIR}/id_rsa" >> /etc/ssh/ssh_config
RUN echo "Host *" >> /etc/ssh/ssh_config && echo " StrictHostKeyChecking no" >> /etc/ssh/ssh_config
RUN chmod -R 600 ${SSHDIR}/* && \
chown -R ${USER}:${USER} ${SSHDIR}/
# check if ssh agent is running or not, if not, run
RUN eval `ssh-agent -s` && ssh-add ${SSHDIR}/id_rsa

EXPOSE 22

#CMD hordesat/hordesat
#CMD supervised-scripts/mpi-run.sh
