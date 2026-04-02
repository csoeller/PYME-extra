"""
CoolLED pE-4000 Serial Driver
================================
Standalone driver + PYME hardware class for the CoolLED pE-4000.

Protocol reference: CoolLED Commands Manual DOC-038
  - Baud rate : 57600
  - Commands  : ASCII, terminated with CR (\r)
  - Responses : terminated with CR+LF (\r\n)
  - Two virtual COM ports appear on Windows; use the LOWER-numbered one
    for commands (the higher one is the event/status port).

Quick usage (standalone):
    dev = PE4000(port='COM3')
    dev.load_wavelength(488)
    dev.set_channel('B', intensity=50, on=True)
    dev.all_off()
    dev.close()
"""

import serial
import time
import logging

log = logging.getLogger(__name__)

# ── channel letters for the pE-4000 ──────────────────────────────────────────
CHANNELS = ('A', 'B', 'C', 'D')

# ── default wavelength layout (can differ by unit config) ────────────────────
DEFAULT_WAVELENGTHS = {
    'A': [365, 385, 405, 435],
    'B': [460, 470, 490, 500],
    'C': [525, 550, 580, 595],
    'D': [635, 660, 740, 770],
}


class PE4000:
    """
    Driver for the CoolLED pE-4000 illumination system.

    Parameters
    ----------
    port : str
        Windows COM port, e.g. 'COM3'.  Use the lower-numbered of the two
        virtual serial ports that appear when the device is connected.
    baud : int
        Baud rate – CoolLED spec is 57600 (default).
    timeout : float
        Read timeout in seconds.
    """

    def __init__(self, port: str, baud: int = 57600, timeout: float = 1.0):
        self._port = port
        self._ser = serial.Serial(
            port=port,
            baudrate=baud,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            timeout=timeout,
        )
        time.sleep(0.1)          # let the port settle
        self._ser.reset_input_buffer()
        log.info("PE4000 connected on %s", port)

        # cache of channel states: {ch: {'on': bool, 'selected': bool, 'intensity': int}}
        self._state = {ch: {'on': False, 'selected': True, 'intensity': 0}
                       for ch in CHANNELS}
        self.refresh_status()    # sync with device

    # ── low-level serial ──────────────────────────────────────────────────────

    def _send(self, cmd: str) -> list[str]:
        """Send a command and collect all response lines."""
        raw = (cmd.strip() + '\r').encode('ascii')
        self._ser.write(raw)
        log.debug("TX: %s", cmd)

        lines = []
        deadline = time.time() + 0.5
        while time.time() < deadline:
            line = self._ser.readline().decode('ascii', errors='replace').strip()
            if line:
                log.debug("RX: %s", line)
                lines.append(line)
                # CoolLED always ends a multi-line response with a CSS… line
                if line.startswith('CSS'):
                    break
        return lines

    # ── status / query ────────────────────────────────────────────────────────

    def refresh_status(self) -> dict:
        """
        Query the device with CSS? and update the internal state cache.
        Returns the parsed state dict.
        """
        lines = self._send('CSS?')
        for line in lines:
            if line.startswith('CSS'):
                self._parse_css(line[3:])   # strip leading 'CSS'
        return self._state

    def _parse_css(self, payload: str):
        """Parse a CSS status payload like AXSF050BSN060... into self._state."""
        i = 0
        while i < len(payload) - 2:
            ch = payload[i].upper()
            if ch not in CHANNELS:
                i += 1
                continue
            sel_char = payload[i + 1].upper()   # S = selected, X = deselected
            on_char  = payload[i + 2].upper()   # N = on,       F = off
            # intensity is the next 2 or 3 digits
            j = i + 3
            while j < len(payload) and payload[j].isdigit():
                j += 1
            intensity = int(payload[i + 3:j]) if j > i + 3 else 0
            self._state[ch] = {
                'selected':  sel_char == 'S',
                'on':        on_char  == 'N',
                'intensity': intensity,
            }
            i = j

    def get_status(self) -> dict:
        """Return the cached channel state (call refresh_status() first if needed)."""
        return dict(self._state)

    def get_version(self) -> list[str]:
        """Return firmware/hardware version strings."""
        return self._send('XVER')

    def get_wavelengths(self) -> dict:
        """
        Query all available wavelengths (LAMBDAS command).
        Returns dict like {'A': [365, 385, 405, 435], 'B': [...], ...}
        """
        lines = self._send('LAMBDAS')
        result = {ch: [] for ch in CHANNELS}
        for line in lines:
            # format: LAMBDA:A0=365
            if line.startswith('LAMBDA:'):
                parts = line[7:]            # 'A0=365'
                ch = parts[0].upper()
                val = parts.split('=')[-1]
                if ch in result:
                    try:
                        result[ch].append(int(val))
                    except ValueError:
                        result[ch].append(val)
        return result

    def get_active_wavelengths(self) -> dict:
        """
        Query wavelengths currently in the ready-for-use position (LAMS).
        Returns dict like {'A': 365, 'B': 460, ...}
        """
        lines = self._send('LAMS')
        result = {}
        for line in lines:
            # format: LAM:A:365
            if line.startswith('LAM:'):
                parts = line.split(':')
                if len(parts) == 3:
                    ch = parts[1].upper()
                    try:
                        result[ch] = int(parts[2])
                    except ValueError:
                        result[ch] = parts[2]
        return result

    # ── wavelength loading ────────────────────────────────────────────────────

    def load_wavelength(self, nm: int) -> list[str]:
        """
        Load a wavelength into the ready-for-use position.
        E.g. load_wavelength(488) → sends 'LOAD:488'
        """
        return self._send(f'LOAD:{nm}')

    # ── channel control ───────────────────────────────────────────────────────

    def set_channel(self, channel: str, intensity: int = None,
                    on: bool = None, selected: bool = True) -> list[str]:
        """
        Set a single channel state.

        Parameters
        ----------
        channel   : 'A', 'B', 'C', or 'D'
        intensity : 0–100 (%)
        on        : True = switch on, False = switch off
        selected  : True = channel participates in global on/off
        """
        ch = channel.upper()
        assert ch in CHANNELS, f"Invalid channel '{ch}'"

        # merge with current cached state
        cur = self._state[ch]
        sel = selected if selected is not None else cur['selected']
        act = on        if on        is not None else cur['on']
        pct = intensity if intensity is not None else cur['intensity']
        pct = max(0, min(100, int(pct)))

        sel_char = 'S' if sel else 'X'
        on_char  = 'N' if act else 'F'
        cmd = f'CSS{ch}{sel_char}{on_char}{pct:03d}'
        lines = self._send(cmd)
        self._state[ch] = {'selected': sel, 'on': act, 'intensity': pct}
        return lines

    def set_intensity(self, channel: str, intensity: int) -> list[str]:
        """Set intensity (0–100 %) for a channel without changing on/off state."""
        return self.set_channel(channel, intensity=intensity)

    def channel_on(self, channel: str) -> list[str]:
        """Turn a single channel on."""
        return self.set_channel(channel, on=True)

    def channel_off(self, channel: str) -> list[str]:
        """Turn a single channel off."""
        return self.set_channel(channel, on=False)

    def all_on(self) -> list[str]:
        """Switch on all *selected* channels (CSN)."""
        lines = self._send('CSN')
        for ch in CHANNELS:
            if self._state[ch]['selected']:
                self._state[ch]['on'] = True
        return lines

    def all_off(self) -> list[str]:
        """Switch off all *selected* channels (CSF)."""
        lines = self._send('CSF')
        for ch in CHANNELS:
            if self._state[ch]['selected']:
                self._state[ch]['on'] = False
        return lines

    def increment_all(self) -> list[str]:
        """Increment intensity of all channels by 1 % (CS+)."""
        return self._send('CS+')

    def decrement_all(self) -> list[str]:
        """Decrement intensity of all channels by 1 % (CS-)."""
        return self._send('CS-')

    # ── control pod lock ──────────────────────────────────────────────────────

    def lock_pod(self) -> list[str]:
        """Disable the physical control pod."""
        return self._send('PORT:P=OFF')

    def unlock_pod(self) -> list[str]:
        """Re-enable the physical control pod."""
        return self._send('PORT:P=ON')

    # ── live reporting ────────────────────────────────────────────────────────

    def live_reporting_on(self) -> list[str]:
        """Enable automatic status reports every ~10 s."""
        return self._send('XLIVE=YES')

    def live_reporting_off(self) -> list[str]:
        """Disable automatic status reports."""
        return self._send('XLIVE=NO')

    # ── cleanup ───────────────────────────────────────────────────────────────

    def close(self):
        """Turn everything off and close the serial port."""
        try:
            self.all_off()
        except Exception:
            pass
        if self._ser.is_open:
            self._ser.close()
        log.info("PE4000 on %s closed", self._port)

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


# =============================================================================
#  PYME hardware interface
# =============================================================================

from PYME.Acquire.Hardware.lasers import Laser

class pe4000laser(Laser):
   def __init__(self, name, turnOn=False, portname='COM3', start_power=0.02):
      self.device = PE4000(portname)
      self.channel = 'C'
      self.name = name
      self.powerControlable = True
      self.device.load_wavelength(560) # check the wavelengths that we actually have
      self.SetPower(start_power) 

   def IsOn(self):
      return self.isOn

   def TurnOn(self):
      self.device.channel_on(self.channel)
      # self.device.set_channel('C', intensity=100.0*self.power, on=True)
      self.isOn = True
      
   def TurnOff(self):
      self.device.channel_off(self.channel)
      # self.device.all_off()
      self.isOn = False

   def SetPower(self, power):
      if power < 0.01 or power > 1:
         raise RuntimeError('Error setting laser power: Power must be between 0 and 1')
      self.power = power
      self.device.set_intensity(self.channel, 100.0*power)

   def GetPower(self):
      return self.power

   def Close(self):
      self.device.channel_off(self.channel)


# =============================================================================
#  Simple CLI / smoke-test
# =============================================================================
if __name__ == '__main__':
    import argparse, sys

    p = argparse.ArgumentParser(description='CoolLED pE-4000 quick test')
    p.add_argument('port', help='COM port, e.g. COM3')
    p.add_argument('--channel', default='A', help='Channel A-D')
    p.add_argument('--intensity', type=int, default=20, help='Intensity 0-100')
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG)

    with PE4000(args.port) as dev:
        print("Version:", dev.get_version())
        print("Active wavelengths:", dev.get_active_wavelengths())
        print(f"\nSetting channel {args.channel} to {args.intensity}%, ON ...")
        dev.set_channel(args.channel, intensity=args.intensity, on=True)
        time.sleep(2)
        print("Turning off ...")
        dev.all_off()
        print("Done.")
