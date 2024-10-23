import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import math
import copy
import yaml
from ..utils.io.yaml_functions import load_yaml

# from .psf.image_psf import load_psf
from ..utils.transform.datatype import *

from ..utils.transform import image_convolution as conv
from .psfs import elliptical_gaussian_3sigmas

from ..utils.io.text import write_txt
from ..utils.io.tiff import write_tif
from ..utils.data_format.visualisation import format_coordinates
from ..utils.visualisation.matplotlib_plots import add_ax_scatter
from ..utils.transform.noise import add_image_noise


class Imager:
    # Optical params
    def __init__(self):
        self.roi_params = {}
        self.roi_params["scale"] = 1e-6
        self.roi_params["reference_point"] = [0, 0]  # center of the imagexy
        self.roi_params["dimension_sizes"] = [1, 1, 1]  # size of the ROI
        self.roi_params["focus_plane"] = 0.0
        self.roi_params["plot_spacing"] = 0.1
        self.roi_params["ranges"] = []
        self.modalities = dict()

        self.field = {}
        self.field["reference_point"] = None
        self.fluorophore_params = (
            dict()
        )  # contiain its photophysical parameters such as excitation and emission
        self.emitters_by_fluorophore = dict()
        self.emitters_by_channel = dict()
        self.writing_dir = ""
        self.identifier = ""

    def set_experiment_name(self, name: str):
        self.identifier = name

    def set_roi_from_file(self, roi_yaml: str):
        with open(roi_yaml, "r") as f:
            roi_p = yaml.safe_load(f)
        self.set_roi_params(**roi_p)

    # methods to get and set roi params
    def set_roi_params(self, **kwargs):
        if kwargs is not None:
            for key, value in kwargs.items():
                self.roi_params[key] = value
            self._calculate_roi_ranges()

    def _calculate_roi_ranges(self):
        dimensions = self.get_roi_params("dimension_sizes")
        reference_xy = self.get_roi_params("reference_point")
        referencez = self.get_roi_params("focus_plane")
        xrange = [
            (reference_xy[0] - (dimensions[0] / 2)),
            (reference_xy[0] + (dimensions[0] / 2)),
        ]
        yrange = [
            (reference_xy[1] - (dimensions[1] / 2)),
            (reference_xy[1] + (dimensions[1] / 2)),
        ]
        zrange = [
            (referencez - (dimensions[2] / 2)),
            (referencez + (dimensions[2] / 2)),
        ]
        self.roi_params["ranges"] = [xrange, yrange, zrange]

    def get_absoulte_reference_point(self):
        focus_plane = self.get_roi_params("focus_plane")
        ref_pt = self.get_roi_params("reference_point")
        ref_pt.append(focus_plane)
        # print(ref_pt)
        return np.array(ref_pt).reshape((1, 3))

    def set_focus(self, focus: float):
        self.roi_params["focus_plane"] = focus
        self._calculate_roi_ranges()

    def set_roi_position(self, x, y):
        self.roi_params["reference_point"] = [x, y]
        self._calculate_roi_ranges()

    def set_roi_sizes(self, x, y):
        self.roi_params["dimension_sizes"][0] = x
        self.roi_params["dimension_sizes"][1] = y
        self._calculate_roi_ranges()
        self._adjust_modalities_imsizes()

    def get_roi_params(self, parameter: str):
        return copy.copy(self.roi_params[parameter])

    def import_field(
        self,
        field_emitters: dict,
        field_scale: float,
        plotting_params: dict,
        reference_point,
        field_sizes,
        **kwargs,
    ):
        if field_scale == self.get_roi_params("scale"):
            print("FOV scale and current scale are the same. No scaling performed")
            scaling_factor = 1
        else:
            scaling_factor = field_scale / self.get_roi_params("scale")
            print(f"scaling factor set to {scaling_factor}")
        self.set_roi_sizes(
            field_sizes[0] * scaling_factor, field_sizes[1] * scaling_factor
        )
        self.field["reference_point"] = reference_point * scaling_factor
        self.recenter_roi()
        for fluoname, emitters in field_emitters.items():
            self.emitters_by_fluorophore[fluoname] = emitters * scaling_factor
        self.plotting_params = plotting_params

    def recenter_roi(self):
        field_xyz = self.field["reference_point"]
        # print(field_xyz)
        self.set_roi_position(field_xyz[0], field_xyz[1])
        self.set_focus(field_xyz[2])

    # fluorophore parameters
    def set_fluorophores_from_file(self, fluo_params: str):
        """
        args
            fluo_params: (string) path to configuration file
        """
        fluo_p = load_yaml(fluo_params)
        self.set_fluorophores_params(**fluo_p)

    def set_fluorophores_params(
        self,
        identifier: str,
        photon_yield: int,
        emission: str,
        blinking_rates: dict,
        **kwargs,
    ):
        fluoname = identifier
        # print(fluo)
        if "plotcolour" in kwargs.keys():
            fluo_color = kwargs["plotcolour"]
        else:
            fluo_color = "blue"
        photons_per_second = photon_yield
        self.fluorophore_params[fluoname] = dict(
            photons_per_second=photons_per_second,
            emission=emission,
            blinking=dict(),
            plotcolour=fluo_color,
        )
        self.set_3state_blinking_params(fluoname, **blinking_rates)

    # self.set_3state_blinking_params()

    def set_fluorophore_photons_per_second(self, fluorophore_name, photon_yield):
        self.fluorophore_params[fluorophore_name]["photons_per_second"] = photon_yield

    def set_3state_blinking_params(
        self,
        fluorophore_name,
        kon,
        koff,
        kbleach,
        initial_state,
        photons_per_blink,
        **kwargs,
    ):
        blink_dictionary = dict(
            kon=kon,
            koff=koff,
            kbleach=kbleach,
            initial_state=initial_state,
            photons_per_blink=photons_per_blink,
        )
        self.fluorophore_params[fluorophore_name]["blinking"] = dict(blink_dictionary)

    # modality module
    def set_imaging_modality_from_file(self, modality_file: str):
        with open(modality_file, "r") as f:
            modality_p = yaml.safe_load(f)
        self.set_imaging_modality(**modality_p)

    def set_imaging_modality(
        self,
        filters: dict,
        psf_params: dict,
        detector: dict,
        emission="blinking",
        modality: str = "modality1",
        **kwargs,
    ):
        self._create_modality(modality)
        self._set_modality_channels(modality, filters)
        self._set_modality_emission(modality, emission)
        self._set_modality_psf(modality, **psf_params)
        self._set_modality_detector(modality, **detector)

    def _create_modality(self, modality: str):
        self.modalities[modality] = dict(
            filters=None, detector=None, psf=None, emission=None
        )

    def _set_modality_channels(self, modality, fluorophores_in_channel):
        self.modalities[modality]["filters"] = fluorophores_in_channel

    def _set_modality_psf(
        self, modality: str, stack_source: str = "generate", **kwargs
    ):
        if stack_source == "generate":
            psf_stack, convolution_type = self._generate_analytical_PSF_stack(**kwargs)
            focus_plane = int((psf_stack.shape)[2] / 2)
            # print(f"focus plane: {focus_plane}")
        else:
            print(f"Loading PSF from file path: {stack_source}")
            psf_stack, convolution_type = self._load_PSF_from_file(
                stack_source, **kwargs
            )
            focus_plane = int(
                (psf_stack.shape)[2] / 2
            )  # this will be changed if the parameter is defined
        # initialising the psf parameters
        self.modalities[modality]["psf"] = dict(
            psf_stack=psf_stack,
            convolution_type=convolution_type,
            focus_plane=focus_plane,
        )
        for key, value in kwargs.items():
            self.modalities[modality]["psf"][key] = value
        # define the psf depth
        if "depth" not in self.modalities[modality]["psf"]:
            depth_default = int((psf_stack.shape)[2] / 2)
            print("No depth parameter found for psf, asigning default")
            self.modalities[modality]["psf"]["depth"] = depth_default

    def _set_modality_emission(self, modality, emission):
        self.modalities[modality]["emission"] = emission

    def _set_modality_detector(
        self,
        modality,
        image_size,  # decpreciated, image size will only be taken from ROI
        pixelsize,
        bits_pixel,
        noise_model=None,
        noise_order=None,
        **kwargs,
    ):
        """
        Set noise and image parameters.
        """
        # image_size_ROI = self.get_roi_params("dimension_sizes")[0:2]
        modality_imsize = self._calculate_imsize_from_ROIranges(pixelsize)
        # print(f"image_size calculated: {modality_imsize} ")
        noise_model_0, noise_order_0 = self._gen_default_noise_model()
        self.modalities[modality]["detector"] = dict(
            image_size=[],
            pixelsize=pixelsize,
            bits_pixel=bits_pixel,
            noise_model=noise_model_0,
            noise_order=noise_order_0,
        )
        if noise_order is not None:
            self.set_noise_order(modality, noise_order)
        if noise_model is not None:
            for noise_t, params in noise_model.items():
                self.set_noise_model_param(modality, noise_t, params)
        self._set_image_size(modality, modality_imsize)

    def _calculate_imsize_from_ROIranges(
        self, mod_pixelsize
    ):  # THere should be specified the scale of the pixelsize of modality
        roi_scale = self.get_roi_params("scale")
        xrange = self.get_roi_params("ranges")[0]
        yrange = self.get_roi_params("ranges")[1]
        xsize = int((xrange[1] - xrange[0]) / mod_pixelsize)
        ysize = int((yrange[1] - yrange[0]) / mod_pixelsize)
        imsize = [xsize, ysize]
        return imsize

    def _adjust_modalities_imsizes(self):
        for mod in self.modalities.keys():
            # print(f"adjusting imsize for {mod}")
            mode_pixelsize = self.modalities[mod]["detector"]["pixelsize"]
            modality_imsize = self._calculate_imsize_from_ROIranges(mode_pixelsize)
            # print(f"image_size calculated: {modality_imsize} ")
            self._set_image_size(mod, modality_imsize)

    def set_noise_order(self, modality, noise_order: list):
        self.modalities[modality]["detector"]["noise_order"] = noise_order
        # print(f"Noise order set as: {noise_order}")

    def set_noise_model_param(self, modality, noise_type, params: dict):
        self.modalities[modality]["detector"]["noise_model"][noise_type] = params
        # print(f"Noise model set as: {noise_type} with params {params}")

    def _gen_default_noise_model(self):
        noise_model = dict(
            binomial={"p": 1},
            gamma={"g": 1},
            gaussian={"sigma": 1},
            poisson={"mock": None},
            conversion={"adu": 1},
            baselevel={"bl": 1},
        )
        noise_order = ["baselevel", "gaussian", "conversion"]
        return noise_model, noise_order

    def _set_image_size(self, modality, imsize):
        # adjust ROI size
        self.modalities[modality]["detector"]["image_size"] = imsize

    # generate and load PSFs
    def _generate_analytical_PSF_stack(
        self, shape=[24, 24, 24], std_devs=[1, 1, 1], **kwargs
    ):
        """
        Generate a eliptical 3D Gaussian PSF
        Shape: in units of pixels
        std_devs: list of standard deviations per dimension
        """
        print(
            f"Generating unitary analytical PSF stack with shape {shape} "
            f"and standard deviations {std_devs}"
        )
        psf_stack = elliptical_gaussian_3sigmas(shape=shape, std_devs=std_devs)
        convolution_type = "3Dvolume"
        return psf_stack, convolution_type

    def _load_PSF_from_file(self, path: str, transpose, unitary, **kwargs):
        # print("loading PSF stack from file")
        psf_stack = load_psf(path, transpose, unitary)
        if unitary:
            convolution_type = "3Dvolume"
        else:
            # convolution_type = "2Dlookup"
            convolution_type = "2Dlookup_optim"
        return psf_stack, convolution_type

    # Imaging methods
    def generate_imaging(
        self,
        modality=None,
        channel="ch0",
        nframes=1,
        nbeads=0,
        save=False,
        noise=False,
        exp_time=1.0,
        **kwargs,
    ):
        """
        Master funciton that generates image sequences depending on the modality
        This function is responsible for
            Preparing the coordinates according to ROI (DONE!)
            Creating the Photons per frame matrix
            Take focus plane, PSF and photon budgets
            Feed these data into the Convolution method (or analytical)
            Adds noise
            Profit!
        """
        if modality is None:
            modality = list(self.modalities.keys())[0]
        # get fluorophores to simulate
        print(f"Simulating imaging from modality: {modality} in channel {channel}")

        fluonames = list(self.modalities[modality]["filters"][channel])
        # prepare a dictionary that conaints the emitters per channel defined
        output_per_fluoname = dict()
        writing_notes = self.identifier + "_" + str(modality) + "_" + str(channel) + "_"
        for fluo in fluonames:  # a channel could capture multiple fluorophores
            # print(fluo)
            writing_notes_fluo = writing_notes + str(fluo)
            emitters = self.get_emitters_in_ROI(fluo)
            n_emitters = emitters.shape[0]
            if n_emitters < 1:
                no_emitters = True
                # if no emitter is in range, an image with noise should be generated
                emitters = np.array(
                    [
                        [0, 0, 0],
                    ]
                )
                photons_frames = np.repeat(0, nframes).reshape((1, nframes))
                emission_notes = "None"
                print(f"photons should be zero: {photons_frames}")
                field_data, psf_data = self._homogenise_scales4convolution_modality(
                    modality, emitters, photons_frames
                )

            else:
                no_emitters = False
                # adjust the absolute XY positions to the ROI dimentions
                emitters[:, 0] = emitters[:, 0] - self.roi_params["ranges"][0][0]
                emitters[:, 1] = emitters[:, 1] - self.roi_params["ranges"][1][0]
                photons_frames, emission_notes = self.calculate_photons_per_frame(
                    modality, fluo, n_emitters, nframes, exp_time
                )
                field_data, psf_data = self._homogenise_scales4convolution_modality(
                    modality, emitters, photons_frames
                )
            # write emitter positions after being placed in the FOV
            gt_notes = writing_notes_fluo + "_usedForImaging"
            emitters_to_export = field_data["field_coordinates"]
            if save:
                self.write_ground_truth_positions(
                    emitters_to_export, "x [nm],y [nm],z [nm]", gt_notes, no_emitters
                )
            convolution_type = self.modalities[modality]["psf"]["convolution_type"]
            psf_size = psf_data["psf_array"].shape
            print(f"size of psf is: {psf_size}")
            simparams = dict(
                field_data=field_data,
                psf_data=psf_data,
                nbeads=nbeads,
                photons_per_bead=self.fluorophore_params[fluo]["blinking"][
                    "photons_per_blink"
                ],
                psf_projection_depth=self.modalities[modality]["psf"]["depth"],
            )
            if convolution_type == "direct":
                pass
            else:
                images, beads = self.images_by_convolutions(
                    convolution_type, **simparams
                )
            writing_notes_fluo = writing_notes_fluo + emission_notes
            # print(f"max and min only photons: {np.max(images)}, {np.min(images)}")
            # Beta implementation to save images before detection
            if save:
                self._save_timeseries_with_beads(images, beads, writing_notes_fluo)

            # wrap up the images from a single fluorophore in the channel
            output_per_fluoname[fluo] = dict(images=images, beads=beads)
        timeseries, beadstack = self._add_fluorophore_signals(output_per_fluoname)
        # # # Up to here only the photon information on arrival
        if noise:
            print("Adding noise")
            timeseries = self._crop_negative(timeseries)
            # print(np.max(timeseries), np.min(timeseries))
            timeseries = self.add_detector_noise(modality, timeseries)
            if beadstack is not None:
                beadstack = self.add_detector_noise(modality, beadstack)
            else:
                beadstack = None
        if save:
            # gt_positions = self.generate_ground_truth_positions(groundtruth_emitters)
            timeseries = self._crop_negative(timeseries)
            beadstack = self._crop_negative(beadstack)
            writing_notes_fluo = writing_notes_fluo + "_withNoise_"

            self._save_timeseries_with_beads(timeseries, beadstack, writing_notes_fluo)
        return timeseries, beadstack

    def write_ground_truth_positions(
        self, emitters, header=None, notes="", no_emitters=False
    ):
        # first, create the list of lines
        textlines = []
        # header = "x [nm],y [nm],z [nm]"
        textlines.append(header)
        if no_emitters:
            self.write_text(textlines, notes)
        else:
            for i in range(emitters.shape[0]):
                string_ = generate_coordinate_textline(
                    emitters[i, 0], emitters[i, 1], emitters[i, 2]
                )
                textlines.append(string_)
            self.write_text(textlines, notes)

    def _crop_negative(self, stack):
        if stack is None:
            return None
        else:
            if np.min(stack) < 0:
                # print("negative numbers found. Setting to zero")
                offset = np.ones(np.shape(stack)) * (-np.min(stack))
                stack = np.add(stack, offset)
            return stack

    def _add_fluorophore_signals(self, output_per_fluoname: dict):
        """
        Dictionary has the fluorophores as keys, each fluorophore
        will have the timeseries and can also have the stack of beads
        if exists
        """
        fluorophores = list(output_per_fluoname.keys())
        shp_ims = output_per_fluoname[fluorophores[0]]["images"].shape
        timeseries = np.zeros((shp_ims))
        for fl in fluorophores:
            both_images = np.array([timeseries, output_per_fluoname[fl]["images"]])
            timeseries = np.sum(both_images, axis=0)
        if output_per_fluoname[fluorophores[0]]["beads"] is not None:
            shp_bds = output_per_fluoname[fluorophores[0]]["beads"].shape
            beadstack = np.zeros((shp_bds))
            for fl in fluorophores:
                both_beads = np.array([beadstack, output_per_fluoname[fl]["beads"]])
                beadstack = np.sum(both_beads, axis=0)
        else:
            beadstack = None
        return timeseries, beadstack

    def add_detector_noise(self, modality, stack):
        """
        INput is an image stack and a dictionary with the
        Types of noise to add, and the parameters needed for each one
        The Noise model also includes an order of noise addition
        """
        # consider all the parameters of the detection system and corrupt
        # the photon data
        noise_params = self.modalities[modality]["detector"]["noise_model"]
        noise_order = self.modalities[modality]["detector"]["noise_order"]
        for noise_type in noise_order:
            # print(noise_type)
            stack = add_image_noise(noise_type, stack, **noise_params[noise_type])
            # includes the processes from collecting the photons, multiplying
            # and all possible source of noise
            # plus the ADU conversion factor
        # Efectively is the digitalisation part and sets the saturation level
        stack = self._adjust_to_pixel_depth(modality, stack)
        return stack

    def _adjust_to_pixel_depth(self, modality, stack):
        bits = self.modalities[modality]["detector"]["bits_pixel"]
        saturaton = (2**bits) - 1
        stack[stack > saturaton] == saturaton
        return stack

    def _save_timeseries_with_beads(self, timeseries_stack, beads_stack, notes):
        notes_images = notes + "_timeseries_photons"
        self.write_tiff_stack(timeseries_stack, notes_images)
        if beads_stack is not None:
            notes_beads = notes + "_calibrationbeads_photons"
            self.write_tiff_stack(beads_stack, notes_beads)

    def images_by_convolutions(
        self,
        convolution_type,
        field_data,
        psf_data,
        nbeads,
        photons_per_bead,
        psf_projection_depth=1,
    ):
        if convolution_type == "2Dlookup":
            # generate convolutions plane by plane then collapse
            downsampled_images, upsampled_images = conv.generate_frames_projection_conv(
                **field_data, **psf_data
            )
            # print("GENERATING BEADS")
            if nbeads >= 1:
                downsampled_beads, upsampled_beads = (
                    conv.generate_beads_projection_conv(
                        **field_data,
                        **psf_data,
                        nbeads=nbeads,
                        bead_photons=photons_per_bead,
                    )
                )
            else:
                downsampled_beads = None
        elif convolution_type == "2Dlookup_optim":
            emitters, images, downsampled_images = (
                conv.generate_frames_projection_conv_optimised2(
                    **field_data, **psf_data
                )
            )
            if nbeads >= 1:
                downsampled_beads, upsampled_beads = (
                    conv.generate_beads_projection_conv(
                        **field_data,
                        **psf_data,
                        nbeads=nbeads,
                        bead_photons=photons_per_bead,
                    )
                )
            else:
                downsampled_beads = None
        elif convolution_type == "3Dvolume":
            # generate convolutions using the whole volume
            downsampled_images = conv.generate_frames_volume_convolution(
                **field_data,
                **psf_data,
                asframes=True,
                psf_projection_depth=psf_projection_depth,
            )
            if nbeads >= 1:
                downsampled_beads = conv.generate_beads_frames_volume_convolution(
                    **field_data,
                    **psf_data,
                    nbeads=nbeads,
                    bead_photons=photons_per_bead,
                    asframes=True,
                )
            else:
                downsampled_beads = None
        return (
            downsampled_images,
            downsampled_beads,
        )  # change at will if needed to get the original images

    def get_emitters_in_ROI(self, fluoname: str):
        # get limits of ROI in xyz
        # to match the previous implementation
        # the limits must be from 0 to maxX in microns
        ranges = self.get_roi_params("ranges")
        points = self._get_emitters_by_fluorophorename(fluoname)
        # print(ranges, points)
        rangesT = np.array(ranges).T
        roi_corners = [rangesT[0].tolist(), rangesT[1].tolist()]
        logical = [inCube(X, roi_corners) for X in points]
        emitters_in_ROI = points[logical, :]
        # up to here the points are only the ones contined
        # in the ROI, but we still need to offset the coordinates
        return copy.copy(emitters_in_ROI)

    def calculate_photons_per_frame(
        self, modality, fluo, n_emitters, nframes, exp_time=1, **kwargs
    ):
        # which model to use
        if n_emitters == 0:
            photon_frames = np.repeat(0, nframes)
            emission_notes = "noemitters"
        else:
            emission = self.modalities[modality]["emission"]
            if emission == "blinking":
                kinetics = dict(self.fluorophore_params[fluo]["blinking"])
                # print(nframes, n_emitters,kinetics)
                photon_frames = self._generate_photons_blinking_modality(
                    nframes, n_emitters, **kinetics
                )
                emission_notes = dictionary2string(kinetics)
            if emission == "constant":
                photons_per_second = self.fluorophore_params[fluo]["photons_per_second"]
                photons_per_frame = exp_time * photons_per_second
                print(f"Average number of photons per frame: {photons_per_frame}")
                photon_frames = self._generate_constant_emission_modality(
                    nframes, n_emitters, photons_per_frame
                )
                emission_notes = "constant_emission_"
        return photon_frames, emission_notes

    def _generate_photons_blinking_modality(self, nframes, nemitters, **kwargs):
        """
        kwargs is intended to hold the kinetic parameters and
        photon budget from the emitter species
        """
        # simulate blinking traces for the N emitters
        # kinetic parameters are in seconds
        tcs = blinking_2states_1bleach_bulk(nframes, nemitters, **kwargs)
        blinking_traces = binary_trace(tcs)  ##
        photons_array = np.random.poisson(
            kwargs["photons_per_blink"], nemitters * nframes
        ).reshape(nemitters, nframes)
        photons_frames = np.multiply(photons_array, blinking_traces)
        return photons_frames

    def _generate_constant_emission_modality(
        self, nframes, nemitters, photons_per_second, exposure_time=1, **kwargs
    ):
        exposure_photons = photons_per_second * exposure_time
        photons_frames = np.random.poisson(
            exposure_photons, nemitters * nframes
        ).reshape(nemitters, nframes)
        return photons_frames

    def _homogenise_scales4convolution_modality(
        self, modality, emitters_in_channel, photons_frames, **kwargs
    ):
        """
        Compiles all data necesary to carry on the generation of frames
        considering the input coordinates
        the calclated photons per frame
        and the PSF of the input modality
        and expresses all dimensions in the dimensions of the PSF

        For example, if the image generator scale is in micrometers
        all the coordinates will be set to nanometers, which is the expected
        scale for the pixelsize of the PSF

        Ideally, the resulting dictionary is passed to a custom convolution
        function

        """
        # check if scaling is needed
        psf_scale = self.modalities[modality]["psf"]["scale"]
        # print(f"psf scale: psf_scale")
        scaling = math.ceil(self.get_roi_params("scale") / psf_scale)
        # print(f"scaling for converting field data for convolution: {scaling}")
        # define field values
        scaling_used = scaling
        field_size = self.modalities[modality]["detector"][
            "image_size"
        ]  # in pixels xyz
        field_pixelsizeXY = (
            self.modalities[modality]["detector"]["pixelsize"] * scaling
        )  # detector pixelsize
        field_coordinates = emitters_in_channel * scaling
        field_zfocus = self.get_roi_params("focus_plane") * scaling
        # ALL PARAMETERS REGARDING FIELD ARE NOW IN PSF SCALE UNITS
        field_scale = psf_scale
        field_data = dict(
            field_scale=field_scale,
            field_size=field_size,
            field_pixelsizeXY=field_pixelsizeXY,
            field_coordinates=field_coordinates,
            field_zfocus=field_zfocus,
            photons_frames=photons_frames,
        )
        # psf parameters
        psf_scale = psf_scale
        psf_array = self.modalities[modality]["psf"]["psf_stack"]
        psf_zstep = self.modalities[modality]["psf"]["voxelsize"][
            2
        ]  # dimensions are assumed to be xyz
        psf_pixelsizeXY = self.modalities[modality]["psf"]["voxelsize"][
            0
        ]  # equal pixelsize laterally
        psf_focus_slice = int(self.modalities[modality]["psf"]["focus_plane"])
        psf_data = dict(
            psf_scale=psf_scale,
            psf_array=psf_array,
            psf_zstep=psf_zstep,
            psf_pixelsizeXY=psf_pixelsizeXY,
            psf_focus_slice=psf_focus_slice,
        )
        # print(f"Matching PSF and coordinates scales")
        return field_data, psf_data

    def _get_emitters_by_fluorophorename(self, fluo):
        return np.array(self.emitters_by_fluorophore[fluo])

    def set_writing_directory(self, saving_dir):
        self.writing_dir = saving_dir

    def write_tiff_stack(self, image_stk, notes):
        path = self.writing_dir
        write_tif(image_stk, path, notes)

    def write_text(self, list_of_lines, notes):
        # text is written line by line
        path = self.writing_dir
        write_txt(list_of_lines, path, notes)

    # methods for visualisation
    def show_field(
        self,
        fluo_type="all",
        view_init=[30, 0, 0],
        initial_pos=True,
        reference_pt=False,
        axesoff=False,
    ):
        # FOCUS
        roi_ranges = self.get_roi_params("ranges")
        z_focus = self.get_roi_params("focus_plane")
        visualisation_scale = 1e-9  # needed for meshgrid
        factor = self.get_roi_params("scale") / visualisation_scale
        xrang = np.array(roi_ranges[0]) * factor
        yrang = np.array(roi_ranges[1]) * factor
        plot_scale = self.get_roi_params("scale") / factor
        # print(f"Showing field on scale: {plot_scale} meters")
        # print(f"dimension ranges in XY pane: {xrang, yrang}")
        # print(f"position of focus {z_focus}")
        xx, yy = np.meshgrid(
            range(int(xrang[0]), int(xrang[1])), range(int(yrang[0]), int(yrang[1]))
        )
        zz = yy * 0 + (z_focus * factor)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(xx, yy, zz, alpha=0.2, cmap="plasma")
        # show ROI reference point
        if reference_pt:
            ref_pt = self.get_absoulte_reference_point() * factor
            add_ax_scatter(ax, format_coordinates(ref_pt))
        # EMITTERS PER FLUOROPHORE SPECIES
        if self.emitters_by_fluorophore is not None:
            for fname, coords in self.emitters_by_fluorophore.items():
                # print(fname)
                add_ax_scatter(
                    ax,
                    format_coordinates(coords * factor, **self.plotting_params[fname]),
                )
        ax.set_box_aspect(
            [ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")]
        )
        ax.view_init(elev=view_init[0], azim=view_init[1], roll=view_init[2])
        if axesoff:
            ax.set_axis_off()
        fig.show
