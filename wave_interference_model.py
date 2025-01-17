import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Define the parameters
rpm = 2500
int_frequency = rpm / 60  # Convert rpm to Hz
omega = 2 * np.pi * int_frequency  # Angular frequency
fint_requency = 3 * rpm / 60  # Frequency of the waves
# Time array
t = np.linspace(0, 0.1, 1000)  # Zoom in by reducing the time range
# Define the waves with harmonics
def generate_wave(t, omega, phase, harmonics=1):
    wave = np.zeros_like(t)
    for n in range(1, harmonics + 1):
        wave+= (1/n) * np.sin(n * omega * t + phase)
        wave += (1/(n+0.5) * np.sin((n+0.5) * omega * t + phase))
    return wave

wave1 = generate_wave(t, omega, 0)
wave2 = generate_wave(t, omega, (1/3)*np.pi)  # Half a harmonic out of phase
wave3 = generate_wave(t, omega, (2/3)*np.pi)  # Quarter a harmonic out of phase

# Define the waves
###wave1 = np.sin(omega * t + (1/3)*np.pi)
#wave2 = np.sin(omega * t + (2/3)*np.pi)  # Half a harmonic out of phase
#wave3 = np.sin(omega * t + (1)*np.pi)  # Quarter a harmonic out of phase

# Sum of the waves
wave_sum = wave1 + wave2 + wave3

# Create the figure and axis
fig, ax = plt.subplots()
line1, = ax.plot(t, wave1, label='Wave 1')
line2, = ax.plot(t, wave2, label='Wave 2')
line3, = ax.plot(t, wave3, label='Wave 3', color='green')  # Third wave in green
line_sum, = ax.plot(t, wave_sum, label='Sum of Waves', color='black')

# Set the axis limits
ax.set_xlim(0, 0.1)  # Zoom in by reducing the x-axis range
ax.set_ylim(-3, 3)

# Add a legend
ax.legend()
# Add the frequency text annotation
frequency_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
# Update function for animation
def update(frame):
    global wave1, wave2, wave3, wave_sum
    frequency = int_frequency * (1 - 0.001 * frame)
    amplitude = 1 - 0.001 * frame
    omega = 2 * np.pi * frequency
    wave1 = amplitude *generate_wave(t, omega,(1/3)*np.pi)
    wave2 = amplitude *generate_wave(t, omega, (2/3)*np.pi)  # Half a harmonic out of phase
    wave3 = amplitude *generate_wave(t, omega, (1)*np.pi)  # Quarter a harmonic out of phase
    #wave1 = amplitude * np.sin(omega * (t + frame / 100) + (1/3) * np.pi)
    #wave2 = amplitude * np.sin(omega * (t + frame / 100) + (2/3) * np.pi)
    #wave3 = amplitude * np.sin(omega * (t + frame / 100) + (2) * np.pi)
    wave_sum = wave1 + wave2 + wave3

    line1.set_ydata(wave1)
    line2.set_ydata(wave2)
    line3.set_ydata(wave3)
    line_sum.set_ydata(wave_sum)

    # Update the frequency text annotation
    frequency_text.set_text(f'Frequency: {frequency:.2f} Hz')

    return line1, line2, line3, line_sum, frequency_text

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=1000, interval=200, blit=True)

# Add the frequency text annotation
frequency_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


# Pause the animation
ani.event_source.stop()

# Display the animation
plt.show()