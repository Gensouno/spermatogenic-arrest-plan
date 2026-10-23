import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ProxyEnvironmentTests(unittest.TestCase):
    def run_config(self, value=None):
        with tempfile.TemporaryDirectory() as folder:
            if value is not None:
                (Path(folder) / '.download-proxy.local').write_text(value)
            env = os.environ.copy()
            env.update(PROJECT_DIR=folder, HTTP_PROXY='http://existing.example:80')
            return subprocess.run(
                ['bash', '-c', 'source "$1"; printf "%s\\n" "$HTTP_PROXY" "${https_proxy:-}" "${NO_PROXY:-}"',
                 'proxy-test', str(ROOT / 'envs/network.sh')],
                env=env, text=True, capture_output=True)

    def test_no_local_config_preserves_environment(self):
        result = self.run_config()
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.splitlines()[0], 'http://existing.example:80')

    def test_local_config_sets_proxy_and_loopback_bypass(self):
        result = self.run_config('http://127.0.0.1:8123\n')
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.splitlines(),
                         ['http://127.0.0.1:8123', 'http://127.0.0.1:8123', 'localhost,127.0.0.1,::1'])

    def test_invalid_scheme_fails_without_direct_fallback(self):
        result = self.run_config('file:///tmp/not-a-proxy')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, '')

    def test_multiline_config_is_not_executed(self):
        result = self.run_config('http://127.0.0.1:8123\nprintf unsafe')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn('unsafe', result.stdout)


if __name__ == '__main__':
    unittest.main()
