from flask import Flask, render_template, request
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
import base64

app = Flask(__name__)

def RungeKutta(t1, t2, h, y0, u0, Torque, J, B1, B2, Ng, Kv, R1, R2, L):
    n = int(np.ceil((t2 - t1) / h)) + 1
    t = np.linspace(t1, t2, n)
    w = np.zeros(n)
    z = np.zeros(n)
    w[0] = y0
    z[0] = u0

    for i in range(n-1):
        k1_w = h * func(t[i], w[i], z[i], Torque, J, B1, B2, Ng, Kv)
        k1_z = h * myfunc2(t[i], z[i], w[i], u0, L, R1, R2, Kv, Ng)

        k2_w = h * func(t[i] + h/2, w[i] + k1_w/2, z[i] + k1_z/2, Torque, J, B1, B2, Ng, Kv)
        k2_z = h * myfunc2(t[i] + h/2, z[i] + k1_z/2, w[i] + k1_w/2, u0, L, R1, R2, Kv, Ng)

        k3_w = h * func(t[i] + h/2, w[i] + k2_w/2, z[i] + k2_z/2, Torque, J, B1, B2, Ng, Kv)
        k3_z = h * myfunc2(t[i] + h/2, z[i] + k2_z/2, w[i] + k2_w/2, u0, L, R1, R2, Kv, Ng)

        k4_w = h * func(t[i] + h, w[i] + k3_w, z[i] + k3_z, Torque, J, B1, B2, Ng, Kv)
        k4_z = h * myfunc2(t[i] + h, z[i] + k3_z, w[i] + k3_w, u0, L, R1, R2, Kv, Ng)

        w[i+1] = w[i] + (k1_w + 2*k2_w + 2*k3_w + k4_w)/6
        z[i+1] = z[i] + (k1_z + 2*k2_z + 2*k3_z + k4_z)/6

    return t, w, z

def func(t, w, z, Torque, J, B1, B2, Ng, Kv):
    Tss = Torque + (0.7 * Torque) * np.sin(np.pi * t)
    return (1 / J) * (Tss - B1 * w - Ng * (z * Kv + B2 * Ng * w))

def myfunc2(t, z, w, u0, L, R1, R2, Kv, Ng):
    return (1 / L) * (-u0 * (R2 + R1) - Kv * Ng * z)

def generate_plot(t_rk, w_rk, z_rk):
  fig, (ax1, ax2) = plt.subplots(2, 1)
  ax1.plot(t_rk, w_rk, '--b', linewidth=2)
  ax1.set_title('Angular Velocity vs. Time')
  ax1.set_xlabel('Time (s)')
  ax1.set_ylabel('Angular Velocity (rad/s)')
  ax1.grid(True)

  ax2.plot(t_rk, z_rk, '--g', linewidth=2)
  ax2.set_title('Current Through The Inductor Over Time')
  ax2.set_xlabel('Time (s)')
  ax2.set_ylabel('Current through Inductor (A)')
  ax2.grid(True)

  # Convert plot to base64 for embedding in HTML
  buffer = BytesIO()
  plt.savefig(buffer, format='png')
  buffer.seek(0)
  plot_data = base64.b64encode(buffer.getvalue()).decode()
  plt.close()

  return plot_data

@app.route('/')
def index():
  return render_template('index.html')

@app.route('/result', methods=['POST'])
def result():
      num_vehicles = int(request.form['num_vehicles'])

      totalWindSpeedProduced = 0
      windSpeed = 4.1667
      airDensity = 1.225
      frontalArea = 2.5
      dragCoefficient = 0.25

      for i in range(num_vehicles):
          vehicleVelocity = float(request.form[f'vehicle_speed_{i+1}'])
          relativeWindSpeed = vehicleVelocity - windSpeed
          dragForce = 0.5 * airDensity * frontalArea * (relativeWindSpeed ** 2) * dragCoefficient
          windSpeedProduced = dragForce / (0.5 * airDensity * frontalArea)
          totalWindSpeedProduced += windSpeedProduced

      print('Total Wind Speed Produced by {} Vehicles: {:.2f} m/s'.format(num_vehicles, totalWindSpeedProduced))
      hours_per_day = 24
      windspeed_mph = totalWindSpeedProduced

      v = windspeed_mph * 0.44704
      Dt = 0.700
      lb = 0.900
      As = Dt * lb
      rho = 1.225
      Pw = 0.5 * rho * As * v**3
      Cp = 0.30
      Pm = Cp * Pw
      
    
      r = Dt / 2
      omega = v / (2 * np.pi * r)
      Torque = Pm / omega
      tl = 0.101
      J = np.pi * (Dt**4 - (Dt - 2 * tl)**4) / 64
      B1 = 0.1
      B2 = 0.1
      N1 = 90
      N2 = 18
      Ng = N1 / N2
      L = 1
      R1 = 8
      R2 = 6
      Kv = 12.5
      t1 = 0
      t2 = 0.01
      y0 = omega
      u0 = 0
      h = 0.0005
      
      # Run the simulation
      # Replace placeholders with appropriate values for your problem
      t_rk, w_rk, z_rk = RungeKutta(t1, t2, h, y0, u0, Torque, J, B1, B2, Ng, Kv, R1, R2, L)
    
      # Generate the plot
      plot_data = generate_plot(t_rk, w_rk, z_rk)
    
      return render_template('result.html', num_vehicles=num_vehicles, Pm=Pm, totalWindSpeedProduced=totalWindSpeedProduced, plot_data=plot_data)


if __name__ == '__main__':
  app.run(host='0.0.0.0')
